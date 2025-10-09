"""
ORF (Open Reading Frame) parser for MetaWorks pipeline

Consolidates functionality from parse_orfs3.py and parse_orfs4.py
"""

from pathlib import Path
from typing import Dict, List, Union, Optional
import sys
import re


class ORFParser:
    """Parse ORFfinder output files"""
    
    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize ORF parser.
        
        Args:
            filepath: Path to ORFfinder output file
        """
        self.filepath = Path(filepath)
        self.orfs: List[Dict] = []
    
    def parse_nucleotide(self) -> List[Dict]:
        """
        Parse nucleotide ORF file.
        
        Returns:
            List of ORF records with nucleotide sequences
        """
        if not self.filepath.exists():
            raise FileNotFoundError(f"ORF file not found: {self.filepath}")
        
        with open(self.filepath, 'r') as f:
            current_orf = None
            
            for line in f:
                line = line.strip()
                
                if line.startswith('>'):
                    if current_orf:
                        self.orfs.append(current_orf)
                    
                    current_orf = self._parse_header(line)
                    current_orf['sequence'] = ''
                
                elif current_orf is not None:
                    current_orf['sequence'] += line
            
            if current_orf:
                self.orfs.append(current_orf)
        
        return self.orfs
    
    def _parse_header(self, header: str) -> Dict:
        """Parse ORF header"""
        id_match = re.match(r'>lcl\|([^_]+)_ORF(\d+)', header)
        if not id_match:
            id_match = re.match(r'>([^\s]+)', header)
            esv_id = id_match.group(1) if id_match else 'unknown'
            orf_num = '1'
        else:
            esv_id = id_match.group(1)
            orf_num = id_match.group(2)
        
        location = re.search(r'\[location=([^\]]+)\]', header)
        length = re.search(r'\[length=(\d+)\]', header)
        
        return {
            'esv_id': esv_id,
            'orf_number': orf_num,
            'location': location.group(1) if location else '',
            'length': int(length.group(1)) if length else 0,
            'full_id': f"{esv_id}_ORF{orf_num}"
        }
    
    def get_longest_orfs(self) -> List[Dict]:
        """Get longest ORF for each ESV"""
        esv_orfs: Dict[str, List[Dict]] = {}
        for orf in self.orfs:
            esv_id = orf['esv_id']
            if esv_id not in esv_orfs:
                esv_orfs[esv_id] = []
            esv_orfs[esv_id].append(orf)
        
        return [max(orf_list, key=lambda x: x['length']) 
                for orf_list in esv_orfs.values()]
    
    def to_fasta(self, orfs: Optional[List[Dict]] = None) -> str:
        """Convert ORFs to FASTA format"""
        if orfs is None:
            orfs = self.orfs
        
        lines = []
        for orf in orfs:
            lines.append(f">{orf['full_id']}")
            lines.append(orf['sequence'])
        
        return '\n'.join(lines)


def parse_orf_output(filepath: Union[str, Path], longest_only: bool = False) -> List[Dict]:
    """Parse ORFfinder output"""
    parser = ORFParser(filepath)
    orfs = parser.parse_nucleotide()
    
    if longest_only:
        orfs = parser.get_longest_orfs()
    
    return orfs


def main():
    """CLI interface"""
    if len(sys.argv) < 2:
        sys.exit("Usage: python -m lib.pseudogene.orf_parser <orf_file> [longest_only]")
    
    filepath = sys.argv[1]
    longest_only = sys.argv[2].lower() == 'true' if len(sys.argv) > 2 else False
    
    try:
        parser = ORFParser(filepath)
        orfs = parser.parse_nucleotide()
        
        if longest_only:
            orfs = parser.get_longest_orfs()
        
        print(parser.to_fasta(orfs))
    except Exception as e:
        sys.exit(f"Error: {e}")


if __name__ == "__main__":
    main()

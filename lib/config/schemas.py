"""Pydantic schemas for configuration validation."""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, validator


class ValidationError(Exception):
    """Raised when configuration validation fails."""

    def __init__(self, errors: List[str]):
        self.errors = errors
        super().__init__(
            f"Validation failed with {len(errors)} error(s):\n"
            + "\n".join(f"  - {e}" for e in errors)
        )


class ConfigError(Exception):
    """Base exception for configuration errors."""

    pass


# ============================================================================
# System Configuration Schema
# ============================================================================


class SystemInfo(BaseModel):
    system: Dict[str, Any]


class SystemPaths(BaseModel):
    paths: Dict[str, str]


class SystemRuntime(BaseModel):
    runtime: Dict[str, Any]


class SystemSchedulers(BaseModel):
    schedulers: Dict[str, Any]


class SystemResources(BaseModel):
    resources: Dict[str, Any]


class SystemLogging(BaseModel):
    logging: Dict[str, Any]


class SystemSecurity(BaseModel):
    security: Dict[str, Any]


class SystemMonitoring(BaseModel):
    monitoring: Dict[str, Any]


class SystemMaintenance(BaseModel):
    maintenance: Dict[str, Any]


class SystemConfig(BaseModel):
    """System configuration schema."""

    system: Dict[str, Any]
    paths: Dict[str, str]
    runtime: Dict[str, Any]
    schedulers: Dict[str, Any]
    resources: Dict[str, Any]
    logging: Dict[str, Any]
    security: Dict[str, Any]
    monitoring: Dict[str, Any]
    maintenance: Dict[str, Any]


# ============================================================================
# Module Configuration Schema
# ============================================================================


class ModuleMetadata(BaseModel):
    name: str
    version: str
    description: str
    author: Optional[str] = None
    enabled_by_default: Optional[bool] = True


class ModuleParameters(BaseModel):
    parameters: Dict[str, Any]


class ModuleFiles(BaseModel):
    files: Dict[str, List[Dict[str, str]]]


class ModuleSoftware(BaseModel):
    software: Dict[str, List[str]]


class ModuleValidation(BaseModel):
    validation: List[Dict[str, Any]]


class ModuleResources(BaseModel):
    resources: Dict[str, Dict[str, int]]


class ModuleDependencies(BaseModel):
    depends_on: List[str] = []


class ModuleCheckpoints(BaseModel):
    checkpoints: List[Dict[str, str]]


class ModuleConfig(BaseModel):
    """Module configuration schema."""

    module: ModuleMetadata
    parameters: ModuleParameters
    files: Optional[ModuleFiles] = None
    software: Optional[ModuleSoftware] = None
    validation: Optional[ModuleValidation] = None
    resources: Optional[ModuleResources] = None
    depends_on: List[str] = []
    checkpoints: Optional[ModuleCheckpoints] = None


# ============================================================================
# User Configuration Schema
# ============================================================================


class PipelineConfig(BaseModel):
    name: str = Field(..., description="Pipeline type (esv or otu)")
    output_dir: str = Field(..., description="Output directory name")
    parallel_jobs: int = Field(default=4, ge=1, le=32)

    @validator("name")
    def validate_pipeline_name(cls, v):
        if v not in ["esv", "otu"]:
            raise ValueError('Pipeline name must be "esv" or "otu"')
        return v

    @validator("output_dir")
    def validate_output_dir(cls, v):
        if not v.replace("_", "").replace("-", "").isalnum():
            raise ValueError(
                "Output directory must contain only alphanumeric characters, underscores, and hyphens"
            )
        return v


class InputConfig(BaseModel):
    sample_source: str = Field(default="folder", description="Sample input method")
    samples_csv: Optional[str] = Field(default=None, description="CSV file with samples")
    fastq_dir: str = Field(..., description="Directory with FASTQ files")

    @validator("sample_source")
    def validate_sample_source(cls, v):
        if v not in ["folder", "csv"]:
            raise ValueError('sample_source must be "folder" or "csv"')
        return v


class ModuleSelection(BaseModel):
    preprocessing: bool = True
    trimming: bool = True
    denoising: bool = True
    classification: bool = True
    pseudogene_filtering: bool = False
    stats: bool = True


class PreprocessingConfig(BaseModel):
    quality_score: int = Field(default=13, ge=0, le=40)
    min_overlap: int = Field(default=25, ge=10, le=100)
    max_mismatch: float = Field(default=0.02, ge=0.0, le=0.5)
    min_match: float = Field(default=0.90, ge=0.0, le=1.0)


class TrimmingConfig(BaseModel):
    adapters: str = Field(..., description="Path to adapter FASTA file")
    min_length: int = Field(default=150, ge=10, le=10000)
    quality_cutoff: str = Field(default="20,20")
    error_rate: float = Field(default=0.1, ge=0.0, le=1.0)
    min_adapter_overlap: int = Field(default=3, ge=1, le=50)
    max_n_bases: int = Field(default=3, ge=0, le=100)
    enable_rc: bool = Field(default=True)


class DenoisingConfig(BaseModel):
    pool_samples: bool = Field(default=True)
    min_cluster_size: int = Field(default=8, ge=1, le=1000)
    threads: int = Field(default=4, ge=1, le=32)


class ClassificationConfig(BaseModel):
    marker: str = Field(default="COI", description="Marker gene")
    memory_gb: int = Field(default=20, ge=1, le=100)
    use_custom_classifier: bool = Field(default=True)
    classifier_path: Optional[str] = Field(default=None)
    builtin_classifier: str = Field(default="fungallsu")
    min_confidence: float = Field(default=0.8, ge=0.0, le=1.0)

    @validator("builtin_classifier")
    def validate_builtin_classifier(cls, v):
        allowed = ["fungallsu", "fungalits_unite", "fungalits_warcup"]
        if v not in allowed:
            raise ValueError(f"builtin_classifier must be one of {allowed}")
        return v


class PseudogeneConfig(BaseModel):
    method: str = Field(default="hmm", description="Filtering method")
    grep_type: int = Field(default=1, ge=1, le=2, description="Grep search type: 1=simple, 2=compound")
    taxon1: str = Field(default="-e Arthropoda", description="First grep pattern")
    taxon2: str = Field(default="", description="Second grep pattern (only for compound search)")
    hmm_profile: str = Field(default="config/hmm/bold.hmm")
    genetic_code: int = Field(default=5, ge=1, le=33)
    orf_start_codon: int = Field(default=2, ge=0, le=2)
    min_orf_length: int = Field(default=30, ge=30, le=10000)
    ignore_nested_orfs: bool = Field(default=True)
    strand: str = Field(default="plus")

    @validator("method")
    def validate_method(cls, v):
        if v not in ["hmm", "orf"]:
            raise ValueError('method must be "hmm" or "orf"')
        return v

    @validator("strand")
    def validate_strand(cls, v):
        if v not in ["both", "plus", "minus"]:
            raise ValueError('strand must be "both", "plus", or "minus"')
        return v


class StatsConfig(BaseModel):
    per_sample: bool = Field(default=True)
    combined: bool = Field(default=True)


class OutputConfig(BaseModel):
    report_type: int = Field(default=1, ge=1, le=2)
    include_intermediate: bool = Field(default=False)
    compress_output: bool = Field(default=True)
    # Note: html_report is reserved for future implementation
    # Currently not implemented, always treated as False


class UserConfig(BaseModel):
    """User configuration schema."""

    pipeline: PipelineConfig
    input: InputConfig
    modules: ModuleSelection
    preprocessing: Optional[PreprocessingConfig] = None
    trimming: Optional[TrimmingConfig] = None
    denoising: Optional[DenoisingConfig] = None
    classification: Optional[ClassificationConfig] = None
    pseudogene_filtering: Optional[PseudogeneConfig] = None
    stats: Optional[StatsConfig] = None
    output: Optional[OutputConfig] = None


# ============================================================================
# Merged Configuration Schema
# ============================================================================


class MergedConfig(BaseModel):
    """Complete merged configuration."""

    system: SystemConfig
    modules: Dict[str, ModuleConfig]
    user: UserConfig

    class Config:
        arbitrary_types_allowed = True

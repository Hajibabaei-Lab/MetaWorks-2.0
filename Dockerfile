# Dockerfile
# docker build -t metaworks .
# docker run -it --rm metaworks bash
# snakemake --jobs 4 --cores 4 --snakefile Snakefile_ESV

FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml \
    && mamba clean --all --yes

RUN wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz \
    && gunzip ORFfinder.gz \
    && chmod +x ORFfinder \
    && mv ORFfinder /opt/conda/envs/MetaWorks/bin/

ARG LD_LIBRARY_PATH
ENV PATH=/opt/conda/envs/MetaWorks/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/conda/envs/MetaWorks/lib:$LD_LIBRARY_PATH

RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /root/.bashrc \
 && echo "conda activate MetaWorks"      >> /root/.bashrc

# Set up working directory
WORKDIR /MetaWorks
COPY . /MetaWorks

CMD ["/bin/bash", "-l"]


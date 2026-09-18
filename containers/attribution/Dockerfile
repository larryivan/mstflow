# Runtime for the source attribution methods that are not packaged in Bioconda:
# FEAST and FastST (R), and the Python libraries SourceID-NMF needs. SourceID-NMF
# itself has no licence, so its code is not copied in here; the workflow fetches it
# from the authors' repository at a pinned commit.
#
#   docker build --platform linux/amd64 -t ghcr.io/larryivan/mstflow-attribution:1.0.0 containers/attribution
FROM rocker/r-ver:4.4.2@sha256:df26749182af64d5263bf64149d51a427b476ed28c4e046997143be3f97fdd7c

# FEAST is CC BY-NC-SA 4.0 (non-commercial); FastST is MIT. Both are fetched at
# fixed commits, and their licence files are kept under /usr/local/share/licenses.
ARG FEAST_COMMIT=2f8f3df8051e0e08341f597a9f4693bfb76b3bf6
ARG FASTST_COMMIT=f59141cee66a07301540da3313c932005f2221d3

LABEL org.opencontainers.image.source="https://github.com/larryivan/mstflow" \
      org.opencontainers.image.description="FEAST, FastST and SourceID-NMF dependencies for MSTFlow"

RUN apt-get update \
    && apt-get install -y --no-install-recommends curl ca-certificates procps python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir --break-system-packages \
        numpy==1.26.4 pandas==2.2.3 scipy==1.14.1 tqdm==4.67.1

# rocker/r-ver installs CRAN packages from a snapshot dated to the R release.
RUN install2.r --error --skipinstalled --ncpus -1 \
        Rcpp RcppArmadillo vegan dplyr reshape2 ggplot2 cowplot ggthemes gridExtra \
        npreg optparse \
    && rm -rf /tmp/downloaded_packages

# The FEAST repository ships compiled objects; remove them so R builds its own.
RUN curl -fsSL "https://github.com/cozygene/FEAST/archive/${FEAST_COMMIT}.tar.gz" | tar -xz -C /tmp \
    && rm -f /tmp/FEAST-${FEAST_COMMIT}/src/*.o /tmp/FEAST-${FEAST_COMMIT}/src/*.so \
    && R CMD INSTALL /tmp/FEAST-${FEAST_COMMIT} \
    && mkdir -p /usr/local/share/licenses/FEAST \
    && cp /tmp/FEAST-${FEAST_COMMIT}/LICENSE.md /usr/local/share/licenses/FEAST/ \
    && rm -rf /tmp/FEAST-${FEAST_COMMIT}

RUN mkdir -p /opt/fastst /usr/local/share/licenses/FastST \
    && curl -fsSL -o /opt/fastst/funST.R \
        "https://raw.githubusercontent.com/joungmin-choi/FastST/${FASTST_COMMIT}/Code/funST.R" \
    && curl -fsSL -o /usr/local/share/licenses/FastST/LICENSE \
        "https://raw.githubusercontent.com/joungmin-choi/FastST/${FASTST_COMMIT}/LICENSE" \
    && printf '%s\n' "${FEAST_COMMIT}" > /usr/local/share/licenses/FEAST/COMMIT \
    && printf '%s\n' "${FASTST_COMMIT}" > /usr/local/share/licenses/FastST/COMMIT

ENV FASTST_HOME=/opt/fastst

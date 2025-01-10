# Env vars

```bash
IMG_NAME=sragent
IMG_VERSION=0.1.1
REGION="us-east1"
GCP_PROJECT_ID="c-tc-429521"
SERVICE_ACCOUNT_EMAIL="nick-nextflow@c-tc-429521.iam.gserviceaccount.com"
SERVICE_ACCOUNT_JSON="c-tc-429521-6f6f5b8ccd93.json"
```

# Docker 

Build the image:

```bash
docker build --platform linux/amd64 \
  -t ${IMG_NAME}:${IMG_VERSION} .
```

Run the image:

```bash
docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v ${PWD}:/data \
  --env DYNACONF="prod" \
  --env EMAIL1="${EMAIL1}" \
  --env EMAIL2="${EMAIL2}" \
  --env NCBI_API_KEY1="${NCBI_API_KEY1}" \
  --env NCBI_API_KEY2="${NCBI_API_KEY2}" \
  --env GCP_SQL_DB_PASSWORD="${GCP_SQL_DB_PASSWORD}" \
  --env OPENAI_API_KEY="${OPENAI_API_KEY}" \
  --platform linux/amd64 \
  ${IMG_NAME}:${IMG_VERSION} --help
```

# GCP Artifact registry

Create, if needed

```bash
DESCRIPTION="SRAgent Docker image"
gcloud artifacts repositories create ${IMG_NAME} \
  --repository-format=docker \
  --project=${GCP_PROJECT_ID} \
  --location=${REGION} \
  --description="${DESCRIPTION}" \
  --async
```

Push

```bash
docker tag ${IMG_NAME}:${IMG_VERSION} \
  ${REGION}-docker.pkg.dev/${GCP_PROJECT_ID}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  && docker push ${REGION}-docker.pkg.dev/${GCP_PROJECT_ID}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION}
```

# Cloud Run jobs

```bash
JOB_NAME="${IMG_NAME}-find-datasets"
gcloud run jobs update ${JOB_NAME} \
  --service-account=${SERVICE_ACCOUNT_EMAIL} \
  --project=${GCP_PROJECT_ID} \
  --region=${REGION} \
  --image=${REGION}-docker.pkg.dev/${GCP_PROJECT_ID}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  --set-env-vars=TZ=America/Los_Angeles \
  --set-env-vars=DYNACONF="prod" \
  --set-env-vars=EMAIL1="nick.youngblut@arcinstitute.org" \
  --set-env-vars=EMAIL2="yusuf.roohani@arcinstitute.org" \
  --set-env-vars=EMAIL3="chris.carpenter@arcinstitute.org" \
  --set-env-vars=EMAIL4="alexander.dobin@arcinstitute.org" \
  --set-env-vars=EMAIL5="hani.goodarzi@arcinstitute.org" \
  --set-secrets=NCBI_API_KEY1=NCBI_API_KEY_NICK:latest \
  --set-secrets=NCBI_API_KEY2=NCBI_API_KEY_YUSUF:latest \
  --set-secrets=NCBI_API_KEY3=NCBI_API_KEY_CHRIS:latest \
  --set-secrets=NCBI_API_KEY4=NCBI_API_KEY_ALEX:latest \
  --set-secrets=NCBI_API_KEY5=NCBI_API_KEY_HANI:latest \
  --set-secrets=GCP_SQL_DB_PASSWORD=GCP_SQL_DB_PASSWORD:latest \
  --set-secrets=OPENAI_API_KEY=OPENAI_API_KEY_SCRECOUNTER:latest \
  --task-timeout=60m \
  --cpu=2 \
  --memory=2Gi \
  --max-retries=0 \
  --args="find-datasets","--no-summaries","Obtain recent single cell RNA-seq datasets in the SRA database"
```

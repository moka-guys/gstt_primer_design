TAG := $(shell git log -1 --pretty=%h)

# Define image names
APP      := gstt_primer_design
REGISTRY := seglh

# Build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(TAG)

.PHONY: push build version cleanbuild

push: build
	docker push $(IMG_VERSIONED)

build: version
	docker buildx build --platform linux/amd64 -t $(IMG_VERSIONED) . || docker build -t $(IMG_VERSIONED) .
	docker tag $(IMG_VERSIONED) $(IMG_VERSIONED) 

cleanbuild:
	docker buildx build --platform linux/amd64 --no-cache -t $(IMG_VERSIONED) .
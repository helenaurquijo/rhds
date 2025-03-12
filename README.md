# rhds

## Setup

To setup this project 

```
cp config-template.env config.env
```

Update the paths in `config.env`.


## Installation

Use `renv` to install relevant packages. In R:

```
install.packages("renv")
renv::restore()
```

## How to run 

### 1. Download the data

```
bash scripts/download-data.sh
Rscript scripts/download-pan-cancer-clinical.r
```


### 2. Extract / clean relevant data

```
Rscript scripts/extract-data.r
```


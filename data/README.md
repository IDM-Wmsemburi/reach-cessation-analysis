# Data Directory Structure

This directory contains the data files required for the REACH cessation analysis. Due to data sensitivity and file size limitations, most raw data files are not included in the repository.

## Directory Structure

```
data/
├── Clean/                     # Processed, analysis-ready data
│   └── data_crt_geospatial.rda   # Main trial dataset 
├── DHS/                       # DHS birth history data 
│   └── *.dta                  # Stata DHS files

├── IHME/                      # IHME vaccination coverage rasters 
│   ├── bcg1.TIF
│   ├── dpt1.TIF
│   ├── dpt3.TIF
│   ├── mcv1.TIF
│   └── polio3.TIF
├── Malaria/                   # Malaria burden data 
│   ├── 2024_GBD2023_Africa_ITN_Use_Rate_2019.TIF
│   ├── 2024_GBD2023_Africa_ITN_Access_2019.TIF
│   ├── 2024_GBD2023_Global_Antimalarial_EFT_2019.TIF
│   ├── 2024_GBD2023_Global_Pf_Incidence_Rate_2019.TIF
│   └── 2024_GBD2023_Global_Pf_Mortality_Rate_2019.TIF
```

## Data Requirements

### 1. REACH Trial Data (`Clean/data_crt_geospatial.rda`)

**Description**: Processed REACH trial consortium data with geospatial coordinates

**Required Variables**:
- `community_id`: Unique cluster identifier
- `trial_phase`: Trial identifier (AVENIR, CHAT, MORDOR I, MORDOR I/II)
- `country`: Country name (Niger, Burkina Faso, Tanzania, Malawi)
- `trt`: Treatment assignment (0 = placebo, 1 = azithromycin)
- `age_start`: Age at start of observation period (months)
- `n`: Number of deaths in observation period
- `person_time`: Total person-time at risk (days)
- `latitude`: Cluster latitude (decimal degrees)
- `longitude`: Cluster longitude (decimal degrees)
- `ihme_u5m_2015`: IHME U5MR estimate for 2015

**Data Sources**:
- AVENIR trial (Niger)
- CHAT trial (Burkina Faso)  
- MORDOR I/II trials (Niger, Tanzania, Malawi)

### 2. DHS Birth History Data (`DHS/*.dta`)

**Description**: Demographic and Health Survey birth history files for hazard estimation

**Expected Files** (not included due to data agreements):
- Burkina Faso DHS files
- Malawi DHS files
- Nigeria DHS files
- Niger DHS files
- Tanzania DHS files

**Required Variables**:
- `v000`: Country code
- `v007`: Survey year
- `v001`: Cluster identifier
- `v002`: Household number
- `v005`: Sample weight
- `v021`: Primary sampling unit
- `v022`: Sample strata
- `v024`: Region
- `v025`: Urban/rural
- Birth history variables (b*, v*)

**Access**: Requires DHS Program data access agreement

### 3. IHME Vaccination Coverage Data (`IHME/*.TIF`)

**Description**: Geospatial rasters of vaccination coverage estimates

**Files**:
- `bcg1.TIF`: BCG first dose coverage
- `dpt1.TIF`: DTP first dose coverage  
- `dpt3.TIF`: DTP third dose coverage
- `mcv1.TIF`: Measles first dose coverage
- `polio3.TIF`: Polio third dose coverage

**Format**: GeoTIFF rasters with global coverage at ~1km resolution

**Source**: Institute for Health Metrics and Evaluation (IHME)

### 4. Malaria Data (`Malaria/*.TIF`)

**Description**: Malaria burden and intervention coverage rasters

**Files**:
- `*_ITN_Use_Rate_2019.TIF`: Insecticide-treated net usage rates
- `*_ITN_Access_2019.TIF`: ITN access rates
- `*_Antimalarial_EFT_2019.TIF`: Antimalarial treatment rates
- `*_Pf_Incidence_Rate_2019.TIF`: P. falciparum incidence per 1,000
- `*_Pf_Mortality_Rate_2019.TIF`: P. falciparum mortality per 100,000

**Format**: GeoTIFF rasters for Africa at ~1km resolution

**Source**: Malaria Atlas Project

## Data Preparation Notes

### Missing Data Handling
- Missing coordinates are excluded from spatial analysis
- Missing person-time observations are excluded
- Missing deaths are treated as zero (after validation)

### Coordinate Systems
- All coordinates should be in WGS84 (EPSG:4326)
- Raster extractions handle coordinate transformations automatically

### Data Quality Checks
The analysis pipeline includes several data validation steps:
- Coordinate bounds checking
- Person-time positivity constraints  
- Treatment coding validation (0/1)
- Age range validation (1-59 months for main analysis)

## Data Access and Ethics

### REACH Trial Data
- Access requires collaboration agreement with trial investigators
- Data use must comply with trial-specific data sharing agreements
- Geographic coordinates may be jittered for privacy

### DHS Data  
- Requires registration and approval from DHS Program
- Subject to DHS data use restrictions
- Geographic data has displacement for confidentiality

### External Covariate Data
- IHME data: Publicly available with citation requirements
- Malaria data: Available through GBD data exchange with registration

## Reproducing the Analysis

To reproduce the analysis with your own data:

1. **Obtain required datasets** following access procedures above
2. **Place files in appropriate directories** following the structure shown
3. **Ensure variable names match** the specifications above
4. **Run data preparation scripts** in the analysis pipeline

For questions about data requirements or access, see the main repository documentation or contact the analysis team.

## File Size Considerations

Due to large file sizes, the following are excluded from version control:
- All raw DHS files (typically 50-200 MB each)
- IHME raster files (typically 100-500 MB each)  
- Malaria raster files (typically 50-300 MB each)
- Large processed datasets (>50 MB)

Consider using Git LFS or external data storage for sharing large datasets.
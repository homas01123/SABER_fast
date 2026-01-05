library(httr2)
library(jsonlite)
library(cli)
library(dplyr)

# ============================================================================
# SENTINEL-2 L1C DOWNLOADER FOR COPERNICUS DATA SPACE
# ============================================================================

# --- SYSTEM CREDENTIALS ---
Sys.setenv(COPERNICUS_USERNAME = "saion523@gmail.com")
Sys.setenv(COPERNICUS_PASSWORD = "SaberHoma$181910")

# ============================================================================
# FUNCTIONS
# ============================================================================

get_token <- function(username, password) {
  cli_h1("Authenticating...")
  
  auth_url <- "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
  
  tryCatch({
    resp <- request(auth_url) %>%
      req_body_form(
        client_id = "cdse-public",
        username = username,
        password = password,
        grant_type = "password"
      ) %>%
      req_perform()
    
    token_data <- resp_body_json(resp)
    cli_alert_success("Authentication successful!")
    return(token_data$access_token)
    
  }, error = function(e) {
    cli_alert_danger("Authentication failed: {e$message}")
    return(NULL)
  })
}

search_products <- function(token, bbox, start_date, end_date, collection, 
                           product_type, max_cloud_cover) {
  cli_h1("Searching for products...")
  
  # Build polygon WKT
  polygon_wkt <- sprintf("POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))",
                        bbox[1], bbox[2],
                        bbox[3], bbox[2],
                        bbox[3], bbox[4],
                        bbox[1], bbox[4],
                        bbox[1], bbox[2])
  
  # Build OData filter
  filter_query <- sprintf(
    "Collection/Name eq '%s' and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and att/OData.CSC.StringAttribute/Value eq '%s') and ContentDate/Start gt %sT00:00:00.000Z and ContentDate/Start lt %sT23:59:59.999Z and OData.CSC.Intersects(area=geography'SRID=4326;%s') and Attributes/OData.CSC.DoubleAttribute/any(att:att/Name eq 'cloudCover' and att/OData.CSC.DoubleAttribute/Value le %f)",
    collection, product_type, start_date, end_date, polygon_wkt, max_cloud_cover
  )
  
  search_url <- "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
  
  tryCatch({
    resp <- request(search_url) %>%
      req_url_query(`$filter` = filter_query, `$top` = 100) %>%
      req_headers(Authorization = paste("Bearer", token)) %>%
      req_perform()
    
    result <- resp_body_json(resp)
    products <- result$value
    
    cli_alert_success("Found {length(products)} products")
    
    if (length(products) > 0) {
      cli_h2("Products:")
      for (i in seq_along(products)) {
        cloud_attr <- Filter(function(x) x$Name == "cloudCover", products[[i]]$Attributes)
        cloud_cover <- if (length(cloud_attr) > 0) round(cloud_attr[[1]]$Value, 2) else "N/A"
        cli_alert_info("{i}. {products[[i]]$Name} - Cloud: {cloud_cover}%")
      }
    }
    
    return(products)
    
  }, error = function(e) {
    cli_alert_danger("Search failed: {e$message}")
    return(list())
  })
}

download_product <- function(token, product, download_dir) {
  product_name <- product$Name
  product_id <- product$Id
  
  cli_h2("Downloading: {product_name}")
  
  if (!dir.exists(download_dir)) {
    dir.create(download_dir, recursive = TRUE)
  }
  
  output_file <- file.path(download_dir, paste0(product_name, ".zip"))
  
  if (file.exists(output_file)) {
    size_mb <- round(file.info(output_file)$size / 1024^2, 2)
    cli_alert_info("Already downloaded ({size_mb} MB)")
    return(TRUE)
  }
  
  download_url <- sprintf(
    "https://zipper.dataspace.copernicus.eu/odata/v1/Products(%s)/$value",
    product_id
  )
  
  tryCatch({
    cli_alert("Starting download (typically 500-800 MB)...")
    cli_alert_info("This may take 5-15 minutes...")
    
    resp <- request(download_url) %>%
      req_headers(Authorization = paste("Bearer", token)) %>%
      req_timeout(3600) %>%
      req_retry(max_tries = 3, is_transient = ~ resp_status(.x) >= 500) %>%
      req_progress() %>%
      req_perform(path = output_file)
    
    size_mb <- round(file.info(output_file)$size / 1024^2, 2)
    cli_alert_success("Complete! ({size_mb} MB)")
    return(TRUE)
    
  }, error = function(e) {
    cli_alert_danger("Download failed: {e$message}")
    if (file.exists(output_file)) file.remove(output_file)
    return(FALSE)
  })
}

download_sentinel2_l1c <- function(username, password, bbox, start_date, end_date,
                                  max_cloud_cover = 10, download_dir = "S2_L1C_Data",
                                  collection = "SENTINEL-2", product_type = "S2MSI1C") {
  
  token <- get_token(username, password)
  if (is.null(token)) return(invisible(NULL))
  
  products <- search_products(token, bbox, start_date, end_date, collection, 
                             product_type, max_cloud_cover)
  
  if (length(products) == 0) {
    cli_alert_warning("No products found")
    return(invisible(NULL))
  }
  
  cli_h1("Downloading {length(products)} products...")
  
  success_count <- 0
  fail_count <- 0
  
  for (i in seq_along(products)) {
    cli_rule(left = sprintf("Product %d/%d", i, length(products)))
    
    success <- download_product(token, products[[i]], download_dir)
    
    if (success) {
      success_count <- success_count + 1
    } else {
      fail_count <- fail_count + 1
    }
    
    if (i < length(products)) Sys.sleep(2)
  }
  
  cli_h1("Summary")
  cli_alert_success("Downloaded: {success_count}/{length(products)}")
  if (fail_count > 0) cli_alert_warning("Failed: {fail_count}/{length(products)}")
  cli_alert_info("Saved to: {download_dir}")
  
  return(invisible(list(success = success_count, failed = fail_count)))
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

USERNAME <- Sys.getenv("COPERNICUS_USERNAME")
PASSWORD <- Sys.getenv("COPERNICUS_PASSWORD")

# --- SEARCH PARAMETERS ---
# Bounding box [minLon, minLat, maxLon, maxLat]
BBOX <- c(-68.532028, 49.01423, -68.106308, 49.193147)

# Date range
START_DATE <- "2019-08-01"
END_DATE <- "2019-08-31"

# Cloud cover maximum (%)
MAX_CLOUD_COVER <- 10

# Collection and product type
COLLECTION <- "SENTINEL-2"
PRODUCT_TYPE <- "S2MSI1C"

# Download directory
DOWNLOAD_DIR <- "S2_L1C_Data"


if (USERNAME == "" || PASSWORD == "") {
  cli_alert_danger("Credentials not set!")
} else {
  results <- download_sentinel2_l1c(
    username = USERNAME,
    password = PASSWORD,
    bbox = BBOX,
    start_date = START_DATE,
    end_date = END_DATE,
    max_cloud_cover = MAX_CLOUD_COVER,
    download_dir = DOWNLOAD_DIR,
    collection = COLLECTION,
    product_type = PRODUCT_TYPE
  )
}

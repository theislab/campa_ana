# Functions to operate on mophological areas identified by pixel-clustering

# convert from dataframe to sparse matrix
to_sparse <- function(d,label_var=NULL) {
  if (!is.null(label_var)) {
    d <- as.matrix(d)
    Matrix::sparseMatrix(i = d[,"y"],j = d[,"x"],x=d[,label_var],index1=FALSE)
  } else {
    Matrix::sparseMatrix(i = d$y,j = d$x,index1=FALSE)
  }
}

to_outlines_sparse <- function(d,size=1) {
  require(imager)
  as.cimg(as.matrix(d)) %>%
    boundary(depth=size) %>%
    where()
}

# locate small regions using EBImage
# TODO: what is the connectivity used here? Seems like 4-connected?
find_small_regions_sparse <- function(d,area=1) {
  label_image <- EBImage::bwlabel(as.matrix(d))
  connected_components <- EBImage::computeFeatures.shape(label_image)
  small_regions <- EBImage::rmObjects(label_image,which(connected_components[,1]>area))
  return(small_regions)
}

find_small_regions_pixset <- function(d,area=1) {
  connected_regions <- as.matrix(d)
  return(small_regions)
}

pad_matrix_zero <- function(m,extra) {
  # add new rows
  new_rows <- matrix(0, ncol=ncol(m), nrow=extra)
  m <- rbind(new_rows, m, new_rows)
  # add new columns
  new_cols <- matrix(0, ncol=extra, nrow=nrow(m))
  m <- cbind(new_cols, m, new_cols)
  return(m)
}


# convert ebimage (matrix) to dataframe using triplet notation of sparse matrix as intermediate
ebimage_to_data_frame <- function(im,labelled=F) {
  triplet <- as(as.matrix(im),"dgTMatrix")
  if(labelled) {
    return(tibble(y=triplet@i,x=triplet@j,id=triplet@x))
  } else {
    return(tibble(y=triplet@i,x=triplet@j))
  }

}

# erode a binary image (provided as a sparse matrix)
erode_sparse <- function(d,width) {
  if (width!=0) {
    return(
      as.matrix(
        EBImage::erode(
          as.matrix(d),
          kern = EBImage::makeBrush(width,shape="diamond")
        )
      )
    )
  } else {
    return(as.matrix(d))
  }
}

# erode a binary image (provided as a sparse matrix)
dilate_sparse <- function(d,width) {
  if (width!=0) {
    # pad matrix before dilation to ensure edges are correct after dilation
    d_pad <- pad_matrix_zero(as.matrix(d),width)
    # dilate
    d_pad_dilated_sparse <- as(
      as.matrix(
        EBImage::dilate(
          d_pad,
          kern = EBImage::makeBrush(width,shape="disc")
        )
      ),"dgTMatrix")
    # adjust coordinates to remove effects of padding
    return(
      as.matrix(
        Matrix::sparseMatrix(i = d_pad_dilated_sparse@i - width,j = d_pad_dilated_sparse@j - width,index1=FALSE,
                             dims=c(max(d_pad_dilated_sparse@i)+1,
                                    max(d_pad_dilated_sparse@j)+1))
      )
    )
  } else {
    return(as.matrix(d))
  }
}


erode_clusters <- function(d,aggregate_by,width=3,remove_small=F) {
  aggregate_by <- enquo(aggregate_by)
  
  # subtract bottom left from coordinates 
  bottom_left <- c(min(d$y),min(d$x))
  
  # nest input using "aggregate_by"
  nested_input <- d %>%
    mutate(y=y-bottom_left[1],x=x-bottom_left[2]) %>%
    select(y,x,!!aggregate_by) %>%
    group_by(!!aggregate_by) %>%
    nest() 
  
  if (remove_small) {
    # locate small regions 
    small_regions <- nested_input %>%
      mutate(sparse = map(data,to_sparse)) %>%
      mutate(filtered = map(sparse,~find_small_regions_sparse(.))) %>%
      mutate(filtered = map(filtered,~ebimage_to_data_frame(.))) %>%
      unnest(filtered) %>%
      select(-sparse,-data) %>%
      ungroup()
    
    # replace small regions with the most common value from neighborhood
    # do not replace pixels where neighborhoods are ambiguous
    replace_values <- small_regions %>%
      select(replace_y=y,replace_x=x) %>%
      mutate(nhd_x = map(replace_x,~seq(.-1,.+1))) %>%
      unnest_longer(nhd_x) %>%
      mutate(nhd_y = map(replace_y,~seq(.-1,.+1))) %>%
      unnest_longer(nhd_y) %>%
      filter(!(replace_y==nhd_y & replace_x==nhd_x)) %>%
      left_join(unnest(nested_input,data),by=c("nhd_x"="x","nhd_y"="y")) %>%
      group_by(replace_y,replace_x) %>%
      add_count(!!aggregate_by) %>%
      filter(n==max(n)) %>%
      distinct(replace_y,replace_x,!!aggregate_by) %>%
      group_by(replace_y,replace_x) %>%
      mutate(n=n()) %>%
      ungroup() %>%
      filter(n==1)
    
    # apply to nested data
    nested_input <- nested_input %>%
      unnest(data) %>%
      # remove small
      anti_join(select(small_regions,y,x), by = c("y", "x")) %>% 
      # replace with non-ambiguous neighbours
      bind_rows(rename(replace_values,y=replace_y,x=replace_x)) %>%
      select(y,x,!!aggregate_by) %>%
      group_by(!!aggregate_by) %>%
      nest() 
  } 
  
  # erode
  res <- nested_input %>%
    mutate(sparse = map(data,to_sparse)) %>%
    mutate(eroded = map(sparse,~erode_sparse(.,width))) %>%
    mutate(eroded = map(eroded,~ebimage_to_data_frame(.))) %>%
    select(-data,-sparse) %>%
    unnest(eroded) %>%
    # add bottom left coordinates back
    mutate(y=y+bottom_left[1],x=x+bottom_left[2])
  
  return(res)
}

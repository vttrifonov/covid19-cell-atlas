#!/bin/bash

s3=s3://$( cat .s3 )
cache=.cache

function s3_exists() {
  local name=$1
  [ ! -z "$( aws s3 ls $s3/$name )" ]
}

function s3_upload() {
  local name=$1
  aws s3 cp - $s3/$name
}

function download1() {
  local name=$1
  local ext=$2
  local url=$( cat input/$name.url )
  [ -f $cache/data/$name.$ext ] && return 0
  mkdir -p $cache/data
  curl -L $url  > $cache/data/$name.$ext
}

download1 mgh h5ad
download1 nih-adaptive h5ad
download1 nih-innate h5ad


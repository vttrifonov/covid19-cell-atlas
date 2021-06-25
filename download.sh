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

function upload1() {
  local name=$1
  local ext=$2
  local url=$( cat input/$name.url )

  s3_exists $name.$ext && return 0
  curl -L $url | s3_upload $name.$ext
}

function download1() {
  local name=$1
  [ -f $cache/$name ] && return 0
  aws s3 cp $s3/$name - > $cache/$name
}

upload1 mgh h5ad
upload1 nih-adaptive h5ad
upload1 nih-innate h5ad

download1 mgh.h5ad
download1 nih-adaptive.h5ad
download1 nih-innate.h5ad


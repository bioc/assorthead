#!/bin/bash

set -e
set -u

mkdir -p sources
mkdir -p licenses
mkdir -p .versions

simple_harvester() {
    local name=$1
    local url=$2
    local version=$3

    local tmpname=sources/${name}
    local vfile=.versions/${name}
    if [ -e ${vfile} ]
    then
        existing_version=$(cat $vfile)
        if [ $existing_version == $version ]
        then
            echo "${name} (${version}) is already present"
            return 0
        fi
    fi

    if [ ! -e $tmpname ]
    then 
        git clone $url $tmpname
    else 
        git -C ${tmpname} fetch --all
    fi

    git -C $tmpname checkout $version
    rm -rf include/$name
    cp -r ${tmpname}/include/$name include/$name

    rm -rf licenses/$name
    mkdir licenses/$name
    cp -r ${tmpname}/LICENSE licenses/${name}

    echo $version > $vfile
}

simple_harvester tatami https://github.com/tatami-inc/tatami v3.0.0
simple_harvester tatami_r https://github.com/tatami-inc/tatami_r v2.0.0
simple_harvester tatami_stats https://github.com/tatami-inc/tatami_stats v1.0.0
simple_harvester manticore https://github.com/tatami-inc/manticore v1.0.2
simple_harvester tatami_chunked https://github.com/tatami-inc/tatami_chunked v2.0.0

simple_harvester byteme https://github.com/LTLA/byteme v1.2.0
simple_harvester aarand https://github.com/LTLA/aarand v1.0.2
simple_harvester powerit https://github.com/LTLA/powerit v2.0.0
simple_harvester WeightedLowess https://github.com/LTLA/CppWeightedLowess v2.0.0
simple_harvester kmeans https://github.com/LTLA/CppKmeans v3.0.1

simple_harvester knncolle https://github.com/knncolle/knncolle v2.0.0
simple_harvester knncolle_annoy https://github.com/knncolle/knncolle_annoy v0.1.0
simple_harvester knncolle_hnsw https://github.com/knncolle/knncolle_hnsw v0.1.0

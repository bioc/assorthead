#!/bin/bash

set -e
set -u

mkdir -p sources
mkdir -p licenses
mkdir -p .versions

get_version_file() {
    echo .versions/$1
}

already_exists() {
    local vfile=$(get_version_file $1)
    if [ -e ${vfile} ]
    then
        existing_version=$(cat $vfile)
        if [ $existing_version == $version ]
        then
            echo 1
            return 0
        fi
    fi

    echo 0
}

git_clone() {
    local name=$1
    local url=$2
    local version=$3

    local tmpname=sources/${name}
    if [ ! -e $tmpname ]
    then 
        git clone $url $tmpname 1>&2
    else 
        git -C ${tmpname} fetch --all 1>&2
    fi

    git -C $tmpname checkout $version 1>&2
    echo $tmpname
}

####################################################

simple_harvester() {
    local name=$1
    local url=$2
    local version=$3

    if [ $(already_exists $name) -eq 1 ]
    then
        echo "${name} (${version}) is already present"
        return 0
    fi

    local tmpname=$(git_clone $name $url $version)
    rm -rf include/$name
    cp -r ${tmpname}/include/$name include/$name

    rm -rf licenses/$name
    mkdir licenses/$name
    cp -r ${tmpname}/LICENSE licenses/${name}

    local vfile=$(get_version_file $name)
    echo $version > $vfile
}

simple_harvester tatami https://github.com/tatami-inc/tatami v3.0.0
simple_harvester tatami_r https://github.com/tatami-inc/tatami_r v2.0.0
simple_harvester tatami_stats https://github.com/tatami-inc/tatami_stats v1.0.0
simple_harvester manticore https://github.com/tatami-inc/manticore v1.0.2
simple_harvester tatami_chunked https://github.com/tatami-inc/tatami_chunked v2.0.0

simple_harvester byteme https://github.com/LTLA/byteme v1.2.2
simple_harvester aarand https://github.com/LTLA/aarand v1.0.2
simple_harvester powerit https://github.com/LTLA/powerit v2.0.0
simple_harvester WeightedLowess https://github.com/LTLA/CppWeightedLowess v2.0.0
simple_harvester kmeans https://github.com/LTLA/CppKmeans v3.0.1

simple_harvester knncolle https://github.com/knncolle/knncolle v2.0.0
simple_harvester knncolle_annoy https://github.com/knncolle/knncolle_annoy v0.1.0
simple_harvester knncolle_hnsw https://github.com/knncolle/knncolle_hnsw v0.1.0

simple_harvester kaori https://github.com/crisprverse/kaori v1.1.2

simple_harvester scran_qc https://github.com/libscran/scran_qc v0.1.0
simple_harvester scran_norm https://github.com/libscran/scran_norm v0.1.0
simple_harvester scran_variances https://github.com/libscran/scran_variances v0.1.0
simple_harvester scran_pca https://github.com/libscran/scran_pca v0.1.0
simple_harvester scran_graph_cluster https://github.com/libscran/scran_graph_cluster v0.1.0
simple_harvester scran_markers https://github.com/libscran/scran_markers v0.1.0
simple_harvester scran_aggregate https://github.com/libscran/scran_aggregate v0.2.0
simple_harvester scran_blocks https://github.com/libscran/scran_blocks v0.1.0

####################################################

annoy_harvester() {
    local name=annoy
    local version=v1.17.2
    local url=https://github.com/spotify/annoy

    if [ $(already_exists $name) -eq 1 ]
    then
        echo "${name} (${version}) is already present"
        return 0
    fi

    local tmpname=$(git_clone $name $url $version)
    rm -rf include/$name
    mkdir include/$name
    cp $tmpname/src/*.h include/$name/

    rm -rf licenses/$name
    mkdir licenses/$name
    cp -r ${tmpname}/LICENSE licenses/${name}

    local vfile=$(get_version_file $name)
    echo $version > $vfile
}

annoy_harvester

####################################################

hnsw_harvester() {
    local name=hnswlib
    local version=v0.8.0
    local url=https://github.com/nmslib/hnswlib

    if [ $(already_exists $name) -eq 1 ]
    then
        echo "${name} (${version}) is already present"
        return 0
    fi

    local tmpname=$(git_clone $name $url $version)
    rm -rf include/$name
    mkdir include/$name
    cp $tmpname/hnswlib/* include/$name/

    rm -rf licenses/$name
    mkdir licenses/$name
    cp -r ${tmpname}/LICENSE licenses/${name}

    local vfile=$(get_version_file $name)
    echo $version > $vfile
}

hnsw_harvester

####################################################

eigen_harvester() {
    local name=Eigen
    local version=3.4.0
    local url=https://gitlab.com/libeigen/eigen

    if [ $(already_exists $name) -eq 1 ]
    then
        echo "${name} (${version}) is already present"
        return 0
    fi

    local tmpname=$(git_clone $name $url $version)
    rm -rf include/$name
    mkdir include/$name
    cp -r $tmpname/Eigen/* include/$name/

    rm -rf licenses/$name
    mkdir licenses/$name
    cp -r ${tmpname}/COPYING.* licenses/${name}

    local vfile=$(get_version_file $name)
    echo $version > $vfile
}

eigen_harvester

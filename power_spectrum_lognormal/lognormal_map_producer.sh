#!/bin/bash

mkdir ./simulations_output/lognormal_map

counter=1
maps_count=1

while [ $counter -le $maps_count ]
do
    flask flask_lognormal.config
    mv ./simulations_output/lognormal_map/map-f1z1.fits ./simulations_output/lognormal_map/lognormal_map.fits 
    mv ./simulations_output/lognormal_map/Xi-f1z1f1z1.dat ./simulations_output/lognormal_map/Xi_flask.dat
    ((counter++))
done

echo $maps_count lognormal map created!

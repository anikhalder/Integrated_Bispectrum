#!/bin/bash

mkdir ../simulations_output/lognormal_maps

counter=1
maps_count=10

while [ $counter -le $maps_count ]
do
    flask flask_lognormal.config
    mv ../simulations_output/lognormal_maps/map-f1z1.fits ../simulations_output/lognormal_maps/lognormal_map_$counter.fits 
    ((counter++))
done

echo $maps_count lognormal maps created!

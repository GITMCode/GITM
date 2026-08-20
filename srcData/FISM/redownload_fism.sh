#!/bin/sh

YEAR=2020
mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat
YEAR=2021
mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat
YEAR=2022
mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat
YEAR=2023
mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat
YEAR=2024
mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat
YEAR=2025
mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat
YEAR=2026
mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat


#YEAR=2017
#mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
#../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
#mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat
#YEAR=2018
#mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
#../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
#mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat
#YEAR=2019
#mv fismflux_daily_${YEAR}.dat fismflux_daily_${YEAR}.save.dat
#../../srcPython/fism2_process.py -start=${YEAR}0101 -end=${YEAR}1231 -euvfile=./euv_59_v2.csv -gitm
#mv fism${YEAR}01_nWaves_059_gitm.dat fismflux_daily_${YEAR}.dat

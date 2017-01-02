#!/bin/bash

less <&0 | cut -d '	' -f 1-3,5,7,9,11,13,15 | tail -n +2 | head -$(($1*($1-1)/2+1)) 

#!/bin/bash

MUTATION_SUMMARY=$1

cut -d ',' -f 3-6 $MUTATION_SUMMARY | tr ',' "\t"

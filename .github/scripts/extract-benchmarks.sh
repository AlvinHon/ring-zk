#!/bin/bash

# Extract benchmark results from Criterion output and format for github-action-benchmark

echo "["

first=true
for dir in target/criterion/*/new/; do
    if [ -f "$dir/estimates.json" ]; then
        benchmark_name=$(basename $(dirname "$dir"))
        mean_ns=$(jq -r '.mean.point_estimate' "$dir/estimates.json")
        
        if [ "$first" = false ]; then
            echo ","
        fi
        first=false
        
        echo "  {"
        echo "    \"name\": \"$benchmark_name\","
        echo "    \"unit\": \"ns\","
        echo "    \"value\": $mean_ns"
        echo -n "  }"
    fi
done

echo ""
echo "]"

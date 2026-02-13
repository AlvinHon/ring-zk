#!/bin/bash
set -euo pipefail

# Extract benchmark results from Criterion output and format for github-action-benchmark
# Note: This script should be committed with executable permissions (chmod +x)

echo "["

first=true
for dir in target/criterion/*/new/; do
    # Skip if glob didn't match any directories
    [ -d "$dir" ] || continue
    
    if [ -f "$dir/estimates.json" ]; then
        benchmark_name=$(basename "$(dirname "$dir")")
        
        # Extract mean value with error handling
        if ! mean_ns=$(jq -r '.mean.point_estimate' "$dir/estimates.json" 2>/dev/null); then
            echo "Error: Failed to parse JSON for $benchmark_name" >&2
            continue
        fi
        
        # Verify we got a valid number
        if [ -z "$mean_ns" ] || ! [[ "$mean_ns" =~ ^[0-9]+\.?[0-9]*$ ]]; then
            echo "Error: Invalid benchmark value for $benchmark_name: $mean_ns" >&2
            continue
        fi
        
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

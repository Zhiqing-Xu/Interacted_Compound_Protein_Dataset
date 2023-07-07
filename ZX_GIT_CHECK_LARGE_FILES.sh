#!/bin/bash

# Set the size limit in bytes (100MB = 104857600 bytes)
sizeLimit=104857600

# Create/overwrite the output file
echo "List of files larger than 100MB" > ZX_GIT_check_large_file_output.txt

# Loop through each file in the current directory and its subdirectories
for file in $(find . -type f); do
    # If the file size is greater than the limit, append its path to the output file
    fileSize=$(stat -c%s "$file")
    if (( fileSize > sizeLimit )); then
        echo $file >> ZX_GIT_check_large_file_output.txt
    fi
done

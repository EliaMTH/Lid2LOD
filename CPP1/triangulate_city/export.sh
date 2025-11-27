#!/bin/bash

basepath="../data/matera"
output_dir="./output"

gui=0

for directory in "${basepath}"/*; do
	if [ -d "$directory" ]; then
		district=$(basename "$directory")
		input="${directory}/output"
		output="${output_dir}/${district}"

		echo "Copying district ${district}: from ${input} to ${output}"
		
		cp -r ${input} ${output}
		echo "Done!"
		echo ""
		echo ""
	fi
done

echo "All Done!"
echo ""


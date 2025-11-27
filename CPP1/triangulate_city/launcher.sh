#!/bin/bash

#exe="build/Qt_6_9_0_for_macOS-Debug/triangulate_city"
exe="build/Desktop_Qt_6_9_0-Debug/triangulate_city"

basepath="../data/matera"
streets="${basepath}/city_street_graph.csv"
gui=0

for directory in "${basepath}"/*; do
	if [ -d "$directory" ]; then
		district=$(basename "$directory")
		ground="${directory}/ground_polygon.off"
		buildings="${directory}/buildings"
		output="${directory}/output"

		echo "Processing district ${district}. Streets: ${streets}, Ground: ${ground}, Buildings: ${buildings}, Output: ${output}"
		./"$exe" "$ground" "$streets" "$buildings" "$output" "$gui"
		echo "Done!"
		echo ""
		echo ""
	fi
done

echo "All Done!"
echo ""

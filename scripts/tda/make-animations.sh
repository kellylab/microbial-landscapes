for x in A B C D E F G
do
		cat subject-trajectories/frames/$x*.png | ffmpeg -framerate 2 -f image2pipe -i - \
				-pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
				subject-trajectories/animations/$x.mp4
done

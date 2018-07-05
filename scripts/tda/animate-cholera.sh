for x in A B C D E F G
do
		ffmpeg -framerate 2 -i cholera-subject-trajectories/frames/$x%d.png \
		-pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
		cholera-subject-trajectories/$x.mp4
done

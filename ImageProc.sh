master="/Users/aesin/Desktop/PUBG_OCR/"
subdir="Lirik/"
video="lirik_clip.mp4"

mkdir -p $master/$subdir/Crop_frames
mkdir -p $master/$subdir/Crop_stacks
mkdir -p $master/$subdir/Block_OCR

width="$(ffprobe -v quiet -print_format json -show_format -show_streams $master/$subdir/$video | grep 'coded_width' | grep -oE '[0-9]+')"
height="$(ffprobe -v quiet -print_format json -show_format -show_streams $master/$subdir/$video | grep 'coded_height' | grep -oE '[0-9]+')"

x_coord=$(echo "$width*0.32" | bc)
y_coord=$(echo "$height*0.67" | bc)
w_size=$(echo "$width*0.35" | bc)
h_size=$(echo "$height*0.055" | bc)


ffmpeg -i $master/$subdir/$video -q:v 2 -vf "crop=$w_size:$h_size:$x_coord:$y_coord:" -r 1 -f image2 /$master/$subdir/Crop_frames/crop_out-%04d.bmp


# List all the .bmp frame files and turn into bash array
frame_files="$(find $master/$subdir/Crop_frames -name "*.bmp" -print0 | xargs -0 ls)"


time {
frame_array=($(echo $frame_files))
# Split these frame files into blocks of 10 which are converted to stacks
block_index=1
block_increment=3

while (( ${#frame_array[@]} != 0 ))
do
	convert ${frame_array[@]::$block_increment} -gravity center -append $master/$subdir/Crop_stacks/Frame_stack_$block_index\.bmp

	convert $master/$subdir/Crop_stacks/Frame_stack_$block_index\.bmp +repage -density 300 -fuzz 90% -negate -resize 350% -define convolve:scale=-100,200%% -morphology Convolve Blur:0x20,100 -level 5% -normalize +contrast -deskew 90% +repage -threshold 15% -depth 1 -bordercolor black -border 1 -fuzz 95% -fill white -draw "color 0,0 floodfill" -alpha off -shave 1x1 -compress Group4 $master/$subdir/Crop_stacks/Frame_stack_$block_index\.tiff

	tesseract -l eng -psm 3 $master/$subdir/Crop_stacks/Frame_stack_$block_index\.tiff $master/$subdir/Block_OCR/Block_OCR_$block_index\.txt


	frame_array=("${frame_array[@]:block_increment}")
	block_index=$[block_index+1]
done
}





overlap_dir="Kill_1"
thresholds=( 15 20 25 30 )

overlap_files="$(find $master/$subdir/$overlap_dir -name "*.bmp" -print0 | xargs -0 ls)"
convert $overlap_files -evaluate-sequence mean $master/$subdir/$overlap_dir/Overlap.bmp

for threshold in "${thresholds[@]}"
do 
	convert $master/$subdir/$overlap_dir/Overlap.bmp +repage -density 300 -fuzz 95% -resize 500% -define convolve:scale=-100,200%% -morphology Convolve Blur:0x30,5 -level 5% -normalize -negate +contrast -deskew 90% +repage -threshold $threshold\% -depth 1 -bordercolor black -border 1 -fuzz 95% -fill white -draw "color 0,0 floodfill" -alpha off -shave 1x1 -compress Group4 $master/$subdir/$overlap_dir/Overlap_$threshold\.tiff

	tesseract -l eng -psm 7 $master/$subdir/$overlap_dir/Overlap_$threshold\.tiff $master/$subdir/$overlap_dir/Overlap_OCR_$threshold\.txt
	
done







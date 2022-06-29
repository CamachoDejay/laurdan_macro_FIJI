/* Simple macro example to split channels and use some simple filtering and segmentation
 * In one of them. This segmented image can later be used for analysis.
 * R. Camacho. [CCI Gothenburg]
 * Current version June 2022, github: camachodejay
*/


close("*");

original_folder = "C:/Users/xcamra/Documents/FIJI_macros/laurdan_macro_FIJI/data/data_3channels/";
split_folder = "C:/Users/xcamra/Documents/FIJI_macros/laurdan_macro_FIJI/data/split_3channels/";

tif_name = "Image 31";

open(original_folder + tif_name + ".tif");
rename(tif_name);
// red is 0
// green is 1
// blue is 2

// getting ch00
ch = "ch00";
ch_name = tif_name + "_" + ch;
selectWindow(tif_name);
run("Duplicate...", "title=["+ ch_name + "] duplicate channels=1");
resetMinAndMax();
saveAs("Tiff", split_folder + ch_name + ".tif");

// getting ch01
ch = "ch01";
ch_name = tif_name + "_" + ch;
selectWindow(tif_name);
run("Duplicate...", "title=["+ ch_name + "] duplicate channels=2");
resetMinAndMax();
saveAs("Tiff", split_folder + ch_name + ".tif");

// getting ch02
ch = "ch02";
ch_name = tif_name + "_" + ch;
selectWindow(tif_name);
run("Duplicate...", "title=["+ ch_name + "] duplicate channels=3");
resetMinAndMax();
saveAs("Tiff", split_folder + ch_name + ".tif");

// getting segmented image
ch = "ch03";
ch_name = tif_name + "_" + ch;

run("Median...", "radius=3");
resetMinAndMax();
setAutoThreshold("Otsu dark");
setOption("BlackBackground", true);
run("Convert to Mask");
saveAs("Tiff", split_folder + ch_name + ".tif");

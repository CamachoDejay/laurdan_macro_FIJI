/* ImageJ macro for GP images analysis
   Quantitative Imaging of Membrane Lipid Order in Cells and Organisms
   Dylan M. Owen, Carles Rentero, Astrid Magenau, Ahmed Abu-Siniyeh and Katharina Gaus
   Nature Protocols 2011
   Original version July 2011
   This macro has been modiffied by R. Camacho. [CCI Gothenburg]
   Current version June 2022, github: camachodejay
*/

/* This macro calculates GP images from Laurdan and di-4-ANEPPDHQ ratiometric images in
   bactch mode (whole chosen folder) obtained using a Leica microscope. The generation
   of HSB images of these GP images has been also implemented.
*/

// house keeping
print("\\Clear");
close("*");
run("Clear Results");
roiManager("reset");
// This makes sure we do not scale images while doing conversions
run("Conversions...", " ");
// to make ife easier I will be handeling general information in a List:
List.clear();
// Keys used in the list 
// DATE and TIME
month_name = "month_name";
day_name = "day_name";
hour_str = "hour_str";
minute_str = "minute_str";
month_id_str = "month_id_str"; // numner of the month formated as string
dayOfMonth_str = "dayOfMonth_str";
year_str = "year_str ";
// USER INPUT
chA = "chA";
chB = "chB";
t = "t";
//lut = "lut";
chC = "chC";
tc = "tc";
Gf = "Gf";
Ques = "Ques";
HSB = "HSB";
lut1 = "lut1";
// folders
results_Dir = "results_Dir";
ordered_images_Dir = "ordered_images_Dir";
disordered_images_Dir = "disordered_images_Dir";
GP_images_Dir = "GP_images_Dir";
histogramGP_Dir = "histogramGP_Dir";
rawGP_images_Dir = "rawGP_images_Dir";
IF_images_Dir = "IF_images_Dir";
GP_IF_images_Dir = "GP_IF_images_Dir";
histogramIF_Dir = "histogramIF_Dir";
HSB_Dir = "HSB_Dir";

// for historical reasons, TODO: check and remove
closeAllImages();

// Select images folder
dir = getDirectory("Choose a Directory ");


// get user input and store it into the list.
getUserInput();

time0 = getTime();
// TODO: setting this to false so I can follow the workflow
setBatchMode(false);

// Here we store in our list the formated strings from date and time
timeDate2List();

// create output directory
// step one get the paths on the List based on date time info
setResultsDirectoriesInList();
// step two we greate the folder that are needed
createOutputFolders();

// Open and save images in 32bit format 
listDir = getFileList(dir);
//var s = 0;
counter = 0;
for (i = 0; i < listDir.length; i++) {
	
	if (endsWith(listDir[i], List.get(chA)+".tif")) {
		open(dir + listDir[i]);
		current_chA = getTitle();
		current_chA = image2gray32bit(current_chA);
		selectWindow(current_chA);
		saveAs("tiff", List.get(ordered_images_Dir) + substring(current_chA,0,lengthOf(current_chA)-9) + "_chA_32bits");
		close();
		counter++; 
		}
		
	if (endsWith(listDir[i], List.get(chB)+".tif")) {
		open(dir + listDir[i]);
		current_chB = getTitle();
		current_chB = image2gray32bit(current_chB);
		selectWindow(current_chB);
		if (List.get(Ques)=="Yes") {
			run("Multiply...","value=" + List.get(Gf));
		}
		saveAs("tiff", List.get(disordered_images_Dir) + substring(current_chB,0,lengthOf(current_chB)-9) + "_chB_32bits");
		close();
		counter++; 
	}
}

if (counter == 0) {
	exit("Error in file inputs, please check that they are tif, and have the proper naming convention");
}
// GP analysis
listOrd = getFileList(List.get(ordered_images_Dir));
listDisord = getFileList(List.get(disordered_images_Dir));

// process disordered ones
for (h = 0, j = h; h < listDisord.length; h++, j++) {
	// Name has the length of # of images
	Names=newArray(listOrd.length);
	
	// open an ordered image, and store the name in the array
	open(List.get(ordered_images_Dir)+listOrd[h]);
	name = getTitle;
	Names[j] = substring(name,0,lengthOf(name)-15);
	
	image_1a = "Image_1a.tif";
	rename(image_1a);
	
	// open a disordered image
	open(List.get(disordered_images_Dir)+listDisord[j]);
	image_2a = "Image_2a.tif"; 
	rename(image_2a);
	
	// calculate difference and sum
	image_diff = "Image_Subs.tif";
	imageCalculator("Substract create 32-bit", image_1a, image_2a);
	rename(image_diff);
	
	image_sum = "Image_Add.tif";
	imageCalculator("Add create 32-bit", image_1a, image_2a);
	rename(image_sum);
	
	// calculate ratio between difference and sum
	GP_image_pre = Names[j] + "_preGP.tif";
	imageCalculator("Divide create 32-bit", image_diff, image_sum);
	setMinAndMax(-1.0000, 1.0000);
	//call("ij.ImagePlus.setDefault16bitRange", 0);
	run(List.get(lut1));
	saveAs("tiff", List.get(rawGP_images_Dir) + GP_image_pre);
	resetMinAndMax();
	
	// Make intensity threshold based on sum image and user defined inpu
	selectImage(image_sum);
	mask_im = "Image_1bit.tif";
	run("Duplicate...", "title=" + mask_im);
	th_value = parseFloat(List.get(t));
	print("the th value is: " + th_value);
	// here the x2 is due to the fact that the mask is defined in the sum image
	setThreshold(th_value*2, 1000000000000000000000000000000.0000);
	setOption("BlackBackground", false);
	run("Convert to Mask");
	// Cahnge to be either 0 or 1
	run("Subtract...", "value=254");
	run("glasbey_on_dark");
	
	GP_image_post = Names[j] + "_GP.tif";
	imageCalculator("Multiply create 32-bit", mask_im, GP_image_pre);
	run(List.get(lut1));
	setMinAndMax(-1.0000, 1.0000);
	saveAs("tiff", List.get(GP_images_Dir) + GP_image_post);
	resetMinAndMax();
	
	// Histogram generation
	histoDir=List.get(histogramGP_Dir);
	HistogramGeneration(GP_image_post, histoDir);
}

// If there is not fluorescent immunostaining
if (List.get(chC) == "none") {

// HSB image generation
	if (List.get(HSB)=="Yes") {
		
		HSBgeneration();
		
	}
// Print information
	listGP = getFileList(List.get(GP_images_Dir));
	time1 = getTime();
	TOTALtime = (time1-time0)/1000;
	printInfo();
}

//TODO: finish this part

// If there is fluorescent immunostaining
else {
// Folder management
	temp_Dir = getDirectory("temp");

	temp_IFmask_Dir = temp_Dir + "IF mask" + File.separator;
	Folder = temp_IFmask_Dir;
	newCleanFolder(Folder);

	File.makeDirectory(List.get(IF_images_Dir));
	File.makeDirectory(List.get(GP_IF_images_Dir));
	File.makeDirectory(List.get(histogramIF_Dir));

// Create a 1-bit mask from Immunofluorescence image and save for each image
	listDir = getFileList(dir);
	Name=newArray(listDir.length);
	counter = 0;
	for (i = 0; i < listDir.length; i++) {
		if (endsWith(listDir[i], chC + ".tif")) {
			open(dir + listDir[i]);
			name = getTitle;
			Name[i]=substring(name,0,lengthOf(name)-9);
			run("32-bit");
			saveAs("tiff",List.get(IF_images_Dir) + Name[i] + "_IF_32bits");
			//run("8-bit");
			setThreshold(List.get(tc),1000000000000000000000000000000.0000);
			run("Convert to Mask");
			run("Divide...","value=255");
			saveAs("tiff", temp_IFmask_Dir + Name[i] + "_GP-IF_1bit");
			close();
			counter++;
		}
	}
	
	if (counter==0) {
		exit("could not find IF images, check naming convention");
	}
	//TODO: finish this part
	exit("here we go");

// GP-IF image
	listTempIF = getFileList(temp_IFmask_Dir);
//	listL = getFileList(GP_images_Dir);
	Name=newArray(listTempIF.length);
	for (h = 0, j = h; h < listTempIF.length; h++, j++) {
		open(temp_IFmask_Dir+listTempIF[h]);
		k = getTitle;
		Name[j]=substring(k,0,lengthOf(k)-15);
		open(GP_images_Dir + Name[j] + "_GP.tif");
		l = getTitle;
		run("8-bit");
		imageCalculator("Multiply create", k, l);
		run(lut1);
		saveAs("tiff", GP_IF_images_Dir + Name[j] + "_GP-IF");

// Histogram and Normal Distribution
		histoDir=histogramIF_Dir;
		HistogramGeneration();
	}

// Folder management
	listDel = getFileList(temp_IFmask_Dir);
	for (i = 0; i < listDel.length; i++) {
		File.delete(temp_IFmask_Dir+listDel[i]); }
	File.delete(temp_IFmask_Dir);

// HSB image generation
	if (HSB=="Yes") {
		HSBgeneration();
	}

// Print information
	listGP = getFileList(List.get(GP_images_Dir));
	time2 = getTime();
	TOTALtime = (time2-time0)/1000;
	printInfo();
}


///////////////// FUNCTIONS ////////////////////



function getUserInput() { 
// get user input and store the values in the List
	// Inicialize choice variables
	CHANNEL1 = newArray("none","ch00","ch01","ch02","ch03","ch04","ch05");
	CHANNEL2 = newArray("ch00","ch01","ch02","ch03","ch04","ch05");
	QUESTION = newArray("Yes","No");
	Lut_Dir = getDirectory("luts");
	lut = getFileList(Lut_Dir);
	
	// Choose image channels and threshold value
	Dialog.create("GP analysis parameters");
	Dialog.addChoice("Acquisition ordered channel:  ", CHANNEL2, "ch00");
	Dialog.addChoice("Acquisition disordered channel:  ", CHANNEL2, "ch01");
	Dialog.addNumber("Lower Threshold Value for GP the mask:  ", 15000);
	Dialog.addChoice("Scale color for GP images:  ", lut, "Rainbow RGB.lut");
	Dialog.addMessage("\n");
	Dialog.addChoice("Immunofluorescence channel:  ", CHANNEL1, "none");
	Dialog.addNumber("Lower Threshold Value for the IF mask:  ", 50);
	Dialog.addMessage("\n");
	Dialog.addNumber("G factor (1 if unknown):  ", 1);
	Dialog.addChoice("Do you want to use G factor for GP image calculation?",QUESTION, "No");
	Dialog.addMessage("\n");
	Dialog.addChoice("Do you want to generate HSB images?",QUESTION, "Yes");
	Dialog.addMessage("\n");
	Dialog.show();
	
	// Feeding variables from dialog choices
	List.set("chA", Dialog.getChoice());
	List.set("chB", Dialog.getChoice());
	List.set("t", Dialog.getNumber());
	//List.set("lut", Dialog.getChoice());
	lut_chosen = Dialog.getChoice();
	List.set("chC", Dialog.getChoice());
	List.set("tc", Dialog.getNumber());
	List.set("Gf", Dialog.getNumber());
	List.set("Ques", Dialog.getChoice());
	List.set("HSB", Dialog.getChoice());
	List.set("lut1",substring(lut_chosen,0,lengthOf(lut_chosen)-4));


}


function timeDate2List() { 
// here we store the values of time date to the list for later use during folder managment
	
	// Folder management
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
	List.set(month_name, MonthNames[month]);
	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
	List.set(day_name, DayNames[dayOfWeek]);
	
	if (hour<10) {hour_tmp = "0"+hour;}
	else {hour_tmp=hour;}
	if (minute<10) {minute_tmp = "0"+minute;}
	else {minute_tmp=minute;}
	if (month<10) {month_id_tmp = "0"+(month+1);}
	else {month_id_tmp=month+1;}
	if (dayOfMonth<10) {dayOfMonth_tmp = "0"+dayOfMonth;}
	else {dayOfMonth_tmp=dayOfMonth;}
	
	List.set(hour_str, hour_tmp);
	List.set(minute_str, minute_tmp);
	List.set(month_id_str, month_id_tmp);
	List.set(dayOfMonth_str, dayOfMonth_tmp);
	List.set(year_str, year);
}

function setResultsDirectoriesInList() { 
// Based on the List info about time and date then we can create all output folder and soter them in the list
	
	tmp_res = dir + "Results "+List.get(year_str)+"-"+List.get(month_id_str)+"-"+List.get(dayOfMonth_str)+" "+List.get(hour_str)+"h"+List.get(minute_str)+ File.separator;
	List.set(results_Dir, tmp_res);
	
	List.set(ordered_images_Dir, List.get(results_Dir) + "Ordered Images" + File.separator);
	
	List.set(disordered_images_Dir, List.get(results_Dir) + "Disordered Images" + File.separator);
	
	List.set(GP_images_Dir, List.get(results_Dir) + "GP images" + File.separator);
	
	List.set(histogramGP_Dir, List.get(GP_images_Dir) + "Histograms" + File.separator);
	
	List.set(rawGP_images_Dir, List.get(results_Dir) + "raw GP images" + File.separator);
	
	List.set(IF_images_Dir, List.get(results_Dir) + "Immunofluorescence Images" + File.separator);
	
	List.set(GP_IF_images_Dir, List.get(results_Dir) + "GP-IF images" + File.separator);
	
	List.set(histogramIF_Dir, List.get(GP_IF_images_Dir) + "Histograms" + File.separator);
	
	List.set(HSB_Dir, List.get(results_Dir) + "HSB images" + File.separator);

}

function createOutputFolders() { 
// once the List is populated by the output folder info then we can create the folders
// it loos silly to create this function but I like to have a clear division of labour
// between geting the paths and then creating the folders
	
	// here we create the output dta structure
	File.makeDirectory(List.get(results_Dir));
	
	newCleanFolder(List.get(ordered_images_Dir));
	
	newCleanFolder(List.get(disordered_images_Dir));
	
	newCleanFolder(List.get(GP_images_Dir));
	
	newCleanFolder(List.get(histogramGP_Dir));
	
	newCleanFolder(List.get(rawGP_images_Dir));

} 

function closeAllImages() {				
// This function closes all images
	while (nImages>0) {
		selectImage(nImages);
		close(); }
}

function newCleanFolder(Folder_path) {					
// This function creates a folder, removing any existing file in a folder with the same name
	File.makeDirectory(Folder_path);
	listFolder = getFileList(Folder_path);
	for (i = 0; i < listFolder.length; i++) {
		File.delete(Folder_path+listFolder[i]); }
}

function image2gray32bit(im_title) {				
// This funcion prepares each image for the analysis
	selectWindow(im_title);
	//run("8-bit");
	run("Grays");
	run("32-bit");
	return im_title;
}

function HistogramGeneration (im_title, histoDir) {
// This funcion obtains the intensity frequency histogram of the images, makes it smoother,
// calculates the normal average distribution and also includes the GP value (and GP value
// corrected by the Gfactor) corresponding to each intensity. An MS Excel file is generated
// with all this data
	selectWindow(im_title);
	tmp_name = substring(im_title, 0, im_title.length - 4);
	
	nBins   = 256;
	Int     = newArray(nBins);
	Cou     = newArray(nBins);
	Smo     = newArray(nBins);
	NAvDist = newArray(nBins);

	getHistogram(values, counts, nBins, -1, 1);
	GP_array  = values;
	GPc_array = newArray(nBins);
	
	if (List.get(Ques)=="Yes") {
		Gf1 = List.get(Gf);
	}else {
		Gf1 = 1;
	}
		
	for (i=0; i<nBins; i++) {
		GP_i = GP_array[i];
		GPc_array[i]=-(1+GP_i+Gf1*GP_i-Gf1)/(-1-GP_i+Gf1*GP_i-Gf1); 
	}
	
	close("*");

	// rolling average over 3
	for (u=0; u<nBins; u++) {
		Int[u]=u;
		Cou[u]=counts[u];
		if (u==0) {
			Smo[u]=counts[u]; 
		}else if (u==255) {
			Smo[u]=counts[u];
		}else {
			Smo[u]=(counts[u-1]+counts[u]+counts[u+1])/3;
		}
	}
	
	
	Array.getStatistics(Cou,min,max,mean,stdDev);
	Sa=(mean*256)-counts[0]-counts[255];
	d=File.open(histoDir + tmp_name+"_Histogram.xls");
	print(d, "Intensity	Counts	Smooth	Norm Av Dist	GP	GP corrected");
	for (k=0; k<256; k++) {
		NAvDist[k]=100*Smo[k]/Sa;
		print(d, Int[k]+"	"+Cou[k]+"	"+Smo[k]+"	"+NAvDist[k]+"	"+GP_array[k]+"	"+GPc_array[k]);
	}
	File.close(d);
}

function HSBgeneration() {				
	// This function generates the HSB image of the GP images as explained in the paper
	closeAllImages();
	setBatchMode(false);

	// Select images folder
	listRAW = getFileList(List.get(rawGP_images_Dir));
	
	BRIGHTNESS = newArray("Order channel","Disorder channel","IF channel");
	BRIGHTNESSnoIF = newArray("Order channel","Disorder channel");
	QUESTION = newArray("Yes","No");

	//TODO: move this path to the list
	//HSB_Dir = List.get(results_Dir) + "HSB images" + File.separator;
	Folder = List.get(HSB_Dir);
	newCleanFolder(Folder);

	Lut_Dir = getDirectory("luts");
	lut = getFileList(Lut_Dir);

	Dialog.create("Choose images and LUT");
	Dialog.addChoice("Select Hue (GP) image: ", listRAW, "none");
	if (List.get(chC) == "none") {
		Dialog.addChoice("Brightness folder: ", BRIGHTNESSnoIF, "Order channel");
	}
	else {
		Dialog.addChoice("Brightness folder: ", BRIGHTNESS, "Order channel");
	}
	Dialog.addMessage("\n");
	Dialog.addChoice("Select color LUT: ",lut, "Rainbow RGB.lut");
	Dialog.addMessage("\n");
	Dialog.addChoice("Do you want to convert the whole folder into HSB images?", QUESTION, "Yes");
	Dialog.show();

// Feeding variables from dialog choices
	H = Dialog.getChoice();
	BRIGHT = Dialog.getChoice();
	Lut = Dialog.getChoice();
	WholeDir=Dialog.getChoice();

	if (BRIGHT=="Order channel") {
		mark = "_chA_32bits.tif";
		brightness_Dir = List.get(results_Dir) + "Ordered Images" + File.separator;
	}
	else if (BRIGHT=="Disorder channel") {
		mark = "_chB_32bits.tif";
		brightness_Dir = List.get(results_Dir) + "Disordered Images" + File.separator;
	}
	else {
		mark = "_IF_32bits.tif";
		brightness_Dir = List.get(results_Dir) + "Immunofluorescence Images" + File.separator;
	}
	
	run("Set Measurements...", "  min limit display redirect=None decimal=5");
	index2 = indexOf(Lut,".lut");
	L = substring(Lut,0,index2);

// HSB whole folder processing
	if (WholeDir == "Yes") {
		for (j = 0; j < listRAW.length; j++) {
			//open(rawGP_images_Dir+listRAW[j]);
			Name = getSingleHSB(List.get(rawGP_images_Dir)+listRAW[j], L, brightness_Dir);
			selectWindow(Name);
			print("where do I save: " + List.get(HSB_Dir) + Name + "_HSB");
			saveAs("tiff", List.get(HSB_Dir) + Name + "_HSB");
			closeAllImages();
		}
	}else {
		Name = getSingleHSB(List.get(rawGP_images_Dir)+H, L, brightness_Dir);
		selectWindow(Name);
		saveAs("tiff", List.get(HSB_Dir) + Name + "_HSB");
		closeAllImages();
	}
	


// Make HSB LUT bar
	HSB_lut = HSB_LUTbar(L);

// This part makes little sense because each image had its own min and max value
/*
// Copy LUT bar to the image
	run("Copy");
	newImage("Panel LUT", "RGB White", 70, 264, 1);
	setBatchMode(false);
	makeRectangle(4,4,32,256);
	run("Paste");
	run("Colors...", "foreground=black background=black selection=yellow");
	run("Line Width...", "line="+2);
	run("Draw");
	run("Select None");
	setFont("Arial", 12);
	setColor(0, 0, 0);
	drawString(d2s(max,2),39,15);
	drawString(d2s(min,2),39,264);
	run("Select None");

	if (WholeDir == "Yes") {
		saveAs("tiff", HSB_Dir + "Lut bar");
		close();
	}
	else {
		saveAs("tiff", HSB_Dir + Name + "_Lut bar");
		open(HSB_Dir + Name + "_HSB.tif");
	}
	*/
	closeAllImages();
}

function HSB_LUTbar(L) { 
// make LUT bar
	setBatchMode(true);
	newImage("Hue", "8-bit Ramp", 256, 256, 1);
	run(L);
	run("Duplicate...", "title=Brightness");
	run("Rotate 90 Degrees Left");
	
	selectWindow("Hue");
	run("RGB Color");
	run("Split Channels");
	
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (red)");
	rename("bR");
	
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (green)");
	rename("bG");
	
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (blue)");
	rename("bB");
	
	run("Merge Channels...", "red=bR green=bG blue=bB gray=*None*");
	rename("HSB LUT");
	selectWindow("Hue (red)");
	close();
	selectWindow("Hue (green)");
	close();
	selectWindow("Hue (blue)");
	close();
	selectWindow("Brightness");
	close();

	selectWindow("HSB LUT");
	run("Rotate 90 Degrees Left");
	run("Size...", "width=32 height=256 interpolation=None");

	setBatchMode(false);
	return "HSB LUT";

}

function getSingleHSB(path2image, L, brightness_Dir) { 
// hadeling single image for HSB
	open(path2image);
	
	name = getTitle();
	Name = substring(name,0,lengthOf(name)-10);
	run(L);
	rename("Hue");
	resetMinAndMax();
	//run("Brightness/Contrast...");
	//waitForUser("set min & max","set min & max");
	getMinAndMax(min,max);
	
	time0 = getTime();

	open(brightness_Dir+Name+mark);
	Bri=getTitle();
	//run("Enhance Contrast", "saturated=0.5 normalize");
	resetMinAndMax();

	selectWindow("Hue");
	run("RGB Color");
	run("Split Channels");

	min_global = 0;
	max_global = 0;
	imageCalculator("Multiply create 32-bit", Bri,"Hue (red)");
	rename("bR");
	resetMinAndMax();
	getMinAndMax(min, max);
	if (min<min_global) {
		min_global = min;
	}
	if (max>max_global) {
		max_global = max;
	}
	//run("8-bit");

	imageCalculator("Multiply create 32-bit", Bri,"Hue (green)");
	rename("bG");
	resetMinAndMax();
	getMinAndMax(min, max);
	if (min<min_global) {
		min_global = min;
	}
	if (max>max_global) {
		max_global = max;
	}
	//run("8-bit")/;

	imageCalculator("Multiply create 32-bit", Bri,"Hue (blue)");
	rename("bB");
	resetMinAndMax();
	getMinAndMax(min, max);
	if (min<min_global) {
		min_global = min;
	}
	if (max>max_global) {
		max_global = max;
	}
	//run("8-bit");

	selectWindow("bR");
	setMinAndMax(min_global, max_global);
	selectWindow("bG");
	setMinAndMax(min_global, max_global);
	selectWindow("bB");
	setMinAndMax(min_global, max_global);
	
	run("Merge Channels...", "c1=bR c2=bG c3=bB create");
    run("RGB Color");
	rename(Name);
	
	return Name;

}

function printInfo() {					// This function prints and saves a summary of the results of the macro
	print("\\Clear");
	print("GP images analysis macro");
	print("   Quantitative Imaging of Membrane Lipid Order in Cells and Organisms");
	print("   Dylan M. Owen, Carles Rentero, Astrid Magenau, Ahmed Abu-Siniyeh and Katharina Gaus");
	print("   Nature Protocols 2011");
	print("   Initial version July 2011");
	print("   Modified by Rafael Camacho, Centre for Cellular Imaging, GU, Sweden");
	print("   Updated on Jun 2022");
	print("\n");
	print("Ordered channel: " + List.get(chA));
	print("Disordered channel: " + List.get(chB));
	print("   Lower threshold value: " + List.get(t));
	if (chC != "none") {
		print("Channel IF: " + List.get(chC));
		print("   Lower threshold value IF mask: " + List.get(tc)); }
	print("");
	print("GP images (" + (listGP.length-1) + ") saved at:");
	print("  " + List.get(GP_images_Dir));
	print("G factor: " + List.get(Gf));
	print("Scale color for GP images: " + List.get(lut1));
	print("");
	if (List.get(chC) != "none") {
		listGP2 = getFileList(GP_IF_images_Dir);
		print("GP-IF images (" + (listGP2.length-1) + ") saved at:");
		print("  " + GP_IF_images_Dir);
		print(""); }
	if (List.get(HSB)=="Yes") {
		print("HSB images saved at:");
		print("  " + List.get(HSB_Dir));
		print("");
	}
	print("");
	print("Execution time: " + d2s(TOTALtime,2) + " seg.");
	print("");
	print("Analysis date: "+List.get(day_name)+", "+List.get(dayOfMonth_str)+" "+List.get(month_name)+" "+List.get(year_str)+" - "+List.get(hour_str)+":"+List.get(minute_str));
	print("ImageJ version " + getVersion());

	selectWindow("Log");
	saveAs("Text", results_Dir + "Log GP images generation "+List.get(year_str)+"-"+List.get(month_id_str)+"-"+List.get(dayOfMonth_str)+" "+List.get(hour_str)+"h"+List.get(minute_str));

	setBatchMode("exit and display");	
}

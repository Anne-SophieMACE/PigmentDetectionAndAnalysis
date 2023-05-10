/*
 * macro detectPigment_quantifySignalAround_vf
 * 
 * This macro requires a multi-channel multi-Z stack, containing at least a transmitted light image. 
 * Pigments will be detected on the trans channel as "minimum" peak spots (= black spots). Pigments 
 * can also be detected in another channel (the red one), if "close enough" to a black pigment they
 * are merged. Then an all other channel expect the nucleus one, an analysis of signal (intensities) 
 * around the pigments (analysis in X,Y and Z) is performed and a results table recapitulates these
 * results.
 * 
 * tested on ImageJ v1.52p on PC & MAC
 * CAREFUL this version requires 1.52p (does not work in version 1.53c or 1.53t because the ROIs are not 
 * correctly associated with the slices!!)
 */

// checks the Fiji version
if( indexOf(getVersion(),"1.52p") == -1 )
	exit("This macro requires Fiji 1.52p");

// closes and resets everything
run("Close All");
print("\\Clear");
roiManager("reset");
roiManager("Show All");
roiManager("Show None");

// asks the user the image he wants to analyze
path = File.openDialog("Choose your composite image"); 
// opens the image
open(path);
//gets image name and directory
img_name = getTitle();
img_path = getDirectory("image");
 // image dimensions: always BF and DAPI at the end, before images to treat: 2 or 3
getDimensions(width, height, channels, slices, frames);

if( frames > slices ){ // if slices and frames are inverted
	run("Properties...", "channels="+channels+" slices="+frames+" frames="=slices);
	Stack.getDimensions(width, height, channels, slices, frames);
}
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); 
run("Set Measurements...", "area mean min centroid integrated redirect=None decimal=5");	

// pauses so that the user checks the slices to keep = the ones not totally out of focus
waitForUser("Look at the image to choose the slices to keep");

// dialog box to ask the microscope used for acquisition
Dialog.create("Acquisition parameters");
Dialog.addNumber("Pixel size (µm)", 0.0645); //0.0647
Dialog.addNumber("Z-step (µm)", 0.2);
Dialog.addNumber("Proportion of positive red pixels:", 0.05);
Dialog.addMessage("Enter the slices you want to study, between 1 and "+slices);
Dialog.addNumber("First slice", 1);
Dialog.addNumber("Last slice", slices);
Dialog.addMessage("Information for pigment detection on trans image:");
Dialog.addNumber("FindMaxima Parameter", 1000);
Dialog.addNumber("Keep black spots with intensity (leave -1 for no threshold) < ", -1);
Dialog.addMessage("Information for pigment detection on red image:");
Dialog.addNumber("FindMaxima Parameter", 1000);
Dialog.addMessage("For the connection between red and trans spots:");
Dialog.addNumber("Distance (radius) in X,Y (µm)", 0.25)
Dialog.addNumber("Distance (radius) in Z (µm)", 0.8);
Dialog.addMessage("For the analysis of proteins around the pigment:");
Dialog.addNumber("Distance (radius) in X,Y (µm)", 0.5);
Dialog.addNumber("Distance (radius) in Z (µm)", 0.4);
Dialog.show();
pix_size_img = Dialog.getNumber();
Zstep = Dialog.getNumber();
prop_red_positive = Dialog.getNumber();
first_slice = Dialog.getNumber();
last_slice = Dialog.getNumber();
prom_param_trans = Dialog.getNumber();
th_black_spot = Dialog.getNumber();
prom_param_red = Dialog.getNumber();
dist_XY_red_blk = Dialog.getNumber();
dist_Z_red_blk = Dialog.getNumber();
dist_XY_prot = Dialog.getNumber();
dist_Z_prot = Dialog.getNumber();

// conversion in pixels/slices
prox_th_XY = 2*dist_XY_prot/pix_size_img;
nbAddZ = floor(dist_Z_prot/Zstep);

// asks the user for the channel correspondances
channel_choices = newArray("Trans","Red","Blue","Green","Cyan (nucleus)");
Dialog.create("Channels correspondance");
for (i = 0; i < channels; i++)
	Dialog.addChoice("Channel "+i+1, channel_choices);
Dialog.addCheckbox("Detect pigments on red channel (if exists)", false);
Dialog.show();

// reads results and assigns the indexes
// -1 if the color does not exist
transPos = -1;
redPos = -1;
greenPos = -1;
bluePos = -1;

for (i_chan = 0; i_chan < channels; i_chan++){
	ans = Dialog.getChoice();
	if( ans == channel_choices[0] )
		transPos = i_chan+1;
	if( ans == channel_choices[1] )
		redPos = i_chan+1;
	if( ans == channel_choices[2] )
		bluePos = i_chan+1;
	if( ans == channel_choices[3] )
		greenPos = i_chan+1;
}

// trans is essential for the analysis
if( transPos == -1 )
	exit("There should be a transmitted light image");

findPigmentRed = Dialog.getCheckbox();

// check if ROIs exist already in the result -> if so increments number for saving
ind_ROI_save = 0;
while( File.exists(img_path+File.separator+"Results_"+img_name+"_"+"_ROI_"+ind_ROI_save+"_"+nbAddZ+"Z_pigment_slices"+first_slice+"-"+last_slice+".xls") )
	ind_ROI_save++;

// duplicates full stack for kept slices
selectWindow(img_name);
run("Duplicate...", "title=img_analysis duplicate slices="+first_slice+"-"+last_slice);
img_name_substack = getTitle();

selectWindow(img_name);
close();

// duplicates trans 
selectWindow(img_name_substack);
run("Duplicate...", "title=trans_img duplicate channels="+transPos);

// ask the user to draw the ROI
setTool("polygon");
waitForUser("Draw the ROI(s) for your study, do not forget to add them in the Manager");

// the user is asked to draw the cells; there shoould be at least 1 ROI
nbCell = roiManager("count");
while(nbCell == 0){
	waitForUser("Draw at least one ROI");
	nbCell = roiManager("count");
}

// renames ROI of cells
for( i_r = 0; i_r < roiManager("count"); i_r++ ){
	roiManager("select", i_r);
	roiManager("rename", "ROI_cell"+i_r+1);
}

// Min projection for trans image
selectWindow("trans_img");
run("Z Project...", "projection=[Min Intensity]");
rename("MinProj");
// finds for each pixel which slice contains the minimum value
findCorrespondanceProjImage("trans_img","MinProj",width, height,slices,"Assoc_slices_trans") ;

// saves the cells as ROI in a temporary file
roiManager("deselect");
roiManager("save",img_path+File.separator+img_name+"_ROI_cell_tmp.zip");

// pigments to detect on the red channel also
if( findPigmentRed && redPos!=-1 ){
	// duplicates the red image
	selectWindow(img_name_substack);
	run("Duplicate...", "title=red_img duplicate channels="+redPos);

	// Max projection for red image
	selectWindow("red_img");
	run("Z Project...", "projection=[Max Intensity]");
	rename("MaxProj");
	// finds for each pixel which slice contains the maximum value
	findCorrespondanceProjImage("red_img","MaxProj",width, height,slices,"Assoc_slices_red") ;
	
	// treats each cell one by one
	for(i_c = 0; i_c < nbCell; i_c ++){
		roiManager("reset");
		
		// keeps only the cell of interest
		roiManager("open", img_path+File.separator+img_name+"_ROI_cell_tmp.zip");
		for( i_r = roiManager("count")-1; i_r >= 0; i_r-- )
			if( i_r != i_c ){
				roiManager("select", i_r);
				roiManager("delete");
			}
		
		// find pigment coordinates in the brightfield: darker points
		coordinates_trans = detectPigmentBrightfield("trans_img","Assoc_slices_trans",prom_param_trans,prox_th_XY,width, height,i_c,th_black_spot);
		
		// coordinates_trans is a table with (X,Y,Z)
		nbPtsTrans = lengthOf(coordinates_trans)/3;
	
		// find pigment coordinates in the red image: lighter points
		selectWindow("red_img");
		roiManager("select", 0);
		detectPigmentRedCh("red_img","Assoc_slices_red",prom_param_red,prox_th_XY,width, height);
		
		spot_red_corres = newArray(nResults);
		
		DrawAndComputeProportionOfSignificantSignal("red_img",prox_th_XY,nbAddZ,prop_red_positive);
	
		// associates red and black spots if they are "close enough"
		dist_same_pig_XY = 2*dist_XY_red_blk/pix_size_img;
		dist_same_pig_Z = dist_Z_red_blk/Zstep;
		
		for(i_sr = 0; i_sr < nResults; i_sr++){
			cx_r = getResult("X", i_sr);
			cy_r = getResult("Y", i_sr);
			cz_r = getResult("Z", i_sr);
			for( i_str = 0; i_str < nbPtsTrans; i_str++){
				if( abs(cx_r-coordinates_trans[i_str*3])<=dist_same_pig_XY && 
					abs(cy_r-coordinates_trans[i_str*3+1])<=dist_same_pig_XY && 
					abs(cz_r-coordinates_trans[i_str*3+2])<=dist_same_pig_Z )
					spot_red_corres[i_sr] = i_str+1;
			}
		}
		
		numberSpotRed = lengthOf(coordinates_trans)/3+1; // increments that starts at the number of spot trans
		spot_red_bool = newArray(nbPtsTrans);
		
		for(i_sr = 0; i_sr < nResults; i_sr ++){
			if( spot_red_corres[i_sr] != 0 ){ // correspondance found
				roiManager("select", spot_red_corres[i_sr]);
				roiManager("rename", "Cell"+i_c+1+"_spot"+spot_red_corres[i_sr]+"_trans-red");
				setResult("Number of spot",i_sr,i_sr+1);
				spot_red_bool[spot_red_corres[i_sr]-1] = 1; // for the final table
			}
			else{ // does not exist: we create a new spot in the Manager 
				selectWindow("red_img");
				Stack.setSlice(getResult("Z",i_sr));
				makeOval(getResult("X",i_sr)-prox_th_XY/2,getResult("Y",i_sr)-prox_th_XY/2, prox_th_XY, prox_th_XY);
				roiManager("add");
				roiManager("select", roiManager("count")-1);
				roiManager("rename", "Cell"+i_c+1+"_spot"+numberSpotRed+"_red");
				numberSpotRed++;
			}
		}
		
		// offers the possibility to remove the red spots 
		waitForUser("Remove the red spots that are not correct");
		numberSpotRed = lengthOf(coordinates_trans)/3+1; // increments that starts at the number of spot trans
		for(i_r = roiManager("count")-1 ; i_r > numberSpotRed; i_r--){
			roiManager("select", i_r);
			roiManager("rename", "Cell"+i_c+1+"_spot"+i_r+"_red");
		}
		// removes the cell to save only pigments for the cell in a temporary file
		roiManager("select",0);
		roiManager("delete");
		roiManager("save",img_path+File.separator+img_name+"_res_cell_"+i_c+"tmp.zip");
	}

	selectWindow("MaxProj");
	close();

	selectWindow("red_img");
	close();
	
	selectWindow("Assoc_slices_red");
	close();
}
else { // no detection of pigment on the red channel
	for(i_c = 0; i_c < nbCell; i_c ++){
		// detects pigment on each cell sequentially
		roiManager("reset");
		roiManager("open", img_path+File.separator+img_name+"_ROI_cell_tmp.zip");
		for( i_r = roiManager("count")-1; i_r >= 0; i_r-- )
			if( i_r != i_c ){ // keeps ONLY the cell of interest
				roiManager("select", i_r);
				roiManager("delete");
			}
		
		// detects pigments on Btightfield
		coordinates_trans = detectPigmentBrightfield("trans_img","Assoc_slices_trans",prom_param_trans,prox_th_XY,width, height,i_c,th_black_spot);
		
		// deletes cell and save temporary zip file with pigments
		roiManager("select", 0);
		roiManager("delete");
		roiManager("save",img_path+File.separator+img_name+"_res_cell_"+i_c+"tmp.zip");
	}
}

selectWindow("MinProj");
close();

selectWindow("trans_img");
close();

selectWindow("Assoc_slices_trans");
close();

// opens (and deletes) temporary ROIs files: cells & pigments for each cell
roiManager("reset");
roiManager("open",img_path+File.separator+img_name+"_ROI_cell_tmp.zip");
del_succ = File.delete(img_path+File.separator+img_name+"_ROI_cell_tmp.zip");
for(i_c = 0; i_c < nbCell; i_c ++){
	roiManager("open",img_path+File.separator+img_name+"_res_cell_"+i_c+"tmp.zip");
	del_succ = File.delete(img_path+File.separator+img_name+"_res_cell_"+i_c+"tmp.zip");
}

// first ROI(s) is/are the ROI(s) drawn by the user
nb_detected_spots = roiManager("count")-nbCell;
run("Clear Results");

selectWindow(img_name_substack);
roiManager("show all");
roiManager("show none");

// green channel exists: we detect signal around pigments in this channel and compute the proportion of
// significant signal (i.e. above the threshold)
if( greenPos != -1 ){
	max_value_green = newArray(nb_detected_spots);
	mean_value_green = newArray(nb_detected_spots);
	max_bg_val_gr = detectSignalAround(img_name_substack,"green_img",greenPos,max_value_green,nb_detected_spots,1/pix_size_img,nbAddZ,mean_value_green,"max",nbCell);
	prop_pos_pix_green = computeProportionOfSignificantSignal("green_img",max_bg_val_gr,nbCell);
}

// red channel exists: we detect signal around pigments in this channel and compute the proportion of
// significant signal (i.e. above the threshold)
if( redPos != -1 ){
	max_value_red = newArray(nb_detected_spots);
	mean_value_red = newArray(nb_detected_spots);
	max_bg_val_red = detectSignalAround(img_name_substack,"red_img",redPos,max_value_red,nb_detected_spots,1/pix_size_img,nbAddZ,mean_value_red,"max",nbCell);
	prop_pos_pix_red = computeProportionOfSignificantSignal("red_img",max_bg_val_red,nbCell);
}

// blue channel exists: we detect signal around pigments in this channel and compute the proportion of
// significant signal (i.e. above the threshold)
if( bluePos != -1 ){
	max_value_blue = newArray(nb_detected_spots);
	mean_value_blue = newArray(nb_detected_spots);
	max_bg_val_bl = detectSignalAround(img_name_substack,"blue_img",bluePos,max_value_blue,nb_detected_spots,1/pix_size_img,nbAddZ,mean_value_blue,"max",nbCell);
	prop_pos_pix_blue = computeProportionOfSignificantSignal("blue_img",max_bg_val_bl,nbCell);
}

// get values on the grey channel: carefull, min to keep, cannot be initialized to 0
min_value_tr = newArray(nb_detected_spots);
for (i = 0; i < nb_detected_spots; i++)
	min_value_tr[i] = 65535;
	
mean_value_tr = newArray(nb_detected_spots);
prop_pos_pix_tr = newArray(nb_detected_spots);
max_bg_val_tr = detectSignalAround(img_name_substack,"trans_img",transPos,min_value_tr,nb_detected_spots,1/pix_size_img,nbAddZ,mean_value_tr,"min",nbCell);

run("Clear Results");
// red channel exists: we save the threshold and close the image
if( redPos != -1 ){
	setResult("Threshold red", 0, max_bg_val_red);
	selectWindow("red_img");
	close();
}

// green channel exists: we save the threshold and close the image
if( greenPos != -1 ){
	setResult("Threshold green", 0, max_bg_val_gr);
	selectWindow("green_img");
	close();
}

// blue channel exists: we save the threshold and close the image
if( bluePos != -1 ){
	setResult("Threshold blue", 0, max_bg_val_bl);
	selectWindow("blue_img");
	close();
}

// results save if at least one channel red, green, blue exists
if( redPos != -1 || greenPos != -1 || bluePos != -1)
	saveAs("Results", img_path+File.separator+"Results_"+img_name+"_thresholds_RGB.xls");
run("Clear Results");

// fill in the final results table
for( i_r = 0; i_r < nb_detected_spots; i_r++ ){
	// spot name as in the Roi Manager
	roiManager("select", i_r+nbCell);
	cur_name = Roi.getName;
	roiManager("select", i_r+nbCell);
	setResult("Spot name",i_r,cur_name);

	// red results
	if( redPos != -1 )
		displayResults("red",i_r,max_value_red[i_r],mean_value_red[i_r],max_bg_val_red,prop_pos_pix_red[i_r]);

	// trans results
	setResult("Min trans around spot",i_r,min_value_tr[i_r]);
	setResult("Mean trans around spot",i_r,mean_value_tr[i_r]);
	if( indexOf(cur_name,"_trans") != -1 )
		setResult("Trans positive",i_r,"yes");
	else 
		setResult("Trans positive",i_r,"no");
	
	if( greenPos != -1 ) // green exists
		displayResults("green",i_r,max_value_green[i_r],mean_value_green[i_r],max_bg_val_gr,prop_pos_pix_green[i_r]);
		
	if( bluePos != -1 ) // blue exists
		displayResults("blue",i_r,max_value_blue[i_r],mean_value_blue[i_r],max_bg_val_bl,prop_pos_pix_blue[i_r]);
	
}
updateResults();

selectWindow("trans_img");
close();

// saves ROI and Results for the image
roiManager("deselect");
roiManager("save",img_path+File.separator+img_name+"_ROI_"+ind_ROI_save+"_"+nbCell+"cells_slices"+first_slice+"-"+last_slice+".zip");
saveAs("Results", img_path+File.separator+"Results_"+img_name+"_"+"_ROI_"+ind_ROI_save+"_"+nbAddZ+"Z_pigment_slices"+first_slice+"-"+last_slice+".xls");


// function to detect in ROIs in the manager if signal is present
// Based on the maximum value (because the protein surrounds the pigment -> cannot take mean), with 
// threshold based on automatic method Otsu
// computes the Otsu threshold (on stack), mean value and max or min (fluo or trans)
function detectSignalAround(img_name,title_img,chan_nb,optimum_value,nb_detected_spots,diam_pix_ana,nbAddZ,mean_value,minOrMax,nbCell){
	selectWindow(img_name);
	run("Duplicate...", "title="+title_img+" duplicate channels="+chan_nb);
	img_chan = getTitle();

	// threshold on the max intensity projection
	run("Z Project...", "projection=[Max Intensity]");
	img_chan_zproj = getTitle(); 
	setAutoThreshold("Triangle dark");
 	getThreshold(max_bg_val,upper);

 	selectWindow(img_chan_zproj);
 	close();

	selectWindow(img_chan);
	resetThreshold();
	selectWindow(img_chan);
	// loop on the spots 
	for( i_r = 0; i_r < nb_detected_spots; i_r++ ){
		run("Clear Results");
		nb_val = 0;
		roiManager("select", i_r+nbCell); // first (nbCell) ROIs are the original ROIs 
		Stack.getPosition(channel, curr_slice, frame);
		run("Measure Stack..."); // measures on the whole stack but keep only on nbAddZ above and below
		for (i_z = maxOf(1,curr_slice-nbAddZ); i_z <= minOf(nSlices,curr_slice+nbAddZ); i_z++){
			if( minOrMax == "max")
				optimum_value[i_r] = maxOf(optimum_value[i_r],getResult("Max", i_z-1));
			if( minOrMax == "min")
				optimum_value[i_r] = minOf(optimum_value[i_r],getResult("Min", i_z-1));
			mean_value[i_r] = mean_value[i_r]+getResult("Mean", i_z-1);
			nb_val++;
		}
		// mean value
		mean_value[i_r] /= nb_val;
	}
	// returns the thresold given by Triangle method
	return max_bg_val;
}
// this function should be called after the FindMaxima on the red channel and computes for each
// peak the proportion of pixels above a threshold (defined with the Triangle method on the max
// projection; coordinates of the points should be in the result table as columns (X,Y,Z)
function DrawAndComputeProportionOfSignificantSignal(img_name,prox_th_XY,nbAddZ,thresh_ratio){
	selectWindow(img_name);
	run("Z Project...", "projection=[Max Intensity]");
	img_chan_zproj = getTitle(); 
	setAutoThreshold("Triangle dark");
 	getThreshold(max_bg_val,upper);

 	selectWindow(img_name);
	roiManager("show all");
	roiManager("show none");
	run("Duplicate...", "title=mask_zstack duplicate");
	setThreshold(max_bg_val,65535);
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Dark black");
	selectWindow("mask_zstack");
	run("Divide...", "value=255 stack"); // division by one so that the sum corresponds to the number of pixels
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

	// reads coordinates of points find by FindMaxima
	tabX_results = newArray(nResults);
	tabY_results = newArray(nResults);
	tabZ_results = newArray(nResults);
	tab_keep = newArray(nResults); // table to keep or not the spot
	for (i = 0; i < nResults; i++) {
		tabX_results[i] = getResult("X", i);
		tabY_results[i] = getResult("Y", i);
		tabZ_results[i] = getResult("Z", i);
	}

	selectWindow("mask_zstack");
	for (i = 0; i < lengthOf(tabX_results); i++) {
		nbPixelPos = 0;
		// circle around the max found for a count using "Measure Stack"
		makeOval(tabX_results[i]-prox_th_XY/2,tabY_results[i]-prox_th_XY/2, prox_th_XY, prox_th_XY);
		roiManager("add");
		roiManager("select", roiManager("count")-1);
		curr_slice = tabZ_results[i];
		run("Clear Results");
		// measures intensity on the slices of interest
		// = current slice + nbAddZ above and below
		run("Measure Stack..."); 
		beg_Z = maxOf(1,curr_slice-nbAddZ);
		end_Z = minOf(nSlices,curr_slice+nbAddZ);
		for (i_z = beg_Z; i_z <= end_Z; i_z++){
			nbPixelPos = nbPixelPos+getResult("RawIntDen", i_z-1);
		}
		nbPixelsROI = getResult("Area", 0)*(end_Z-beg_Z+1);
		
		// keeps the red pigment only if a certain number of values are
		// above the threshold 
		if( nbPixelPos/nbPixelsROI > thresh_ratio )
			tab_keep[i] = 1;
		
		// removes the circle 
		roiManager("select", roiManager("count")-1);
		roiManager("delete");
	}
	
	// writes in the result table the coordinates of the kept points
	run("Clear Results");
	ind_res = 0;
	for (i = 0; i < lengthOf(tabX_results); i++) {
		if( tab_keep[i] == 1 ) {
			setResult("X", ind_res, tabX_results[i]);
			setResult("Y", ind_res, tabY_results[i]);
			setResult("Z", ind_res, tabZ_results[i]);
			ind_res++;
		}
	}
	selectWindow("mask_zstack");
	close();
}

// this function computes the proportion of significant signal (= above the threshold given by threshold_img)
// on the image names img_name
function computeProportionOfSignificantSignal(img_name,threshold_img,nbCell){
	// mask image, threshold given by the previous function
	selectWindow(img_name);
	roiManager("show all");
	roiManager("show none");
	// creates mask for all values above threshold_img
	run("Duplicate...", "title=mask_zstack duplicate");
	setThreshold(threshold_img,65535);
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Dark black");
	selectWindow("mask_zstack");
	run("Divide...", "value=255 stack"); // division by one so that the sum corresponds to the number of pixels
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // scale in pixels so that the area is in pixel

	prop_table = newArray(roiManager("count"));
	selectWindow("mask_zstack");
	// loop on spots of the cell
	for( i_r = 0; i_r < nb_detected_spots; i_r++ ){
		nbPixelPos = 0;
		roiManager("Select", i_r+nbCell);
		Stack.getPosition(channel, curr_slice, frame);
		run("Clear Results");
		// measure on the slices of interest
		run("Measure Stack..."); 
		beg_Z = maxOf(1,curr_slice-nbAddZ);
		end_Z = minOf(nSlices,curr_slice+nbAddZ);
		for (i_z = beg_Z; i_z <= end_Z; i_z++) // counts number of positive pixels on the ROI
			nbPixelPos = nbPixelPos+getResult("RawIntDen", i_z-1);
		// counts total number of pixels
		nbPixelsROI = getResult("Area", 0)*(end_Z-beg_Z+1);
		prop_table[i_r] = nbPixelPos/nbPixelsROI;
	}

	selectWindow("mask_zstack");
	close();

	return prop_table;
}

// function to detect the pigments based on FindMaxima function on the Min projection
// the user can remove manually some spots 
function detectPigmentBrightfield(tit_img,tit_assoc_slice,prom_param_trans,prox_th_XY,width,height,num_cell,th_blk){

	selectWindow("MinProj");
	roiManager("select", 0);
	run("Find Maxima...", "prominence="+prom_param_trans+" light output=List");

	val_z = newArray(nResults);
	
	selectWindow(tit_assoc_slice);
	for(i_r = 0; i_r < nResults; i_r++)
		val_z[i_r] = getPixel(getResult("X",i_r),getResult("Y",i_r));
	
	// wait because sometimes image was not selected!!
	selectWindow(tit_img); 
	wait(10);
	// creates round ROIs around each spot found by FindMaxima
	for(i_r = 0; i_r < nResults; i_r++){
		Stack.setSlice(val_z[i_r]);
		makeOval(getResult("X",i_r)-prox_th_XY/2,getResult("Y",i_r)-prox_th_XY/2, prox_th_XY, prox_th_XY);
		roiManager("add");
		roiManager("select", roiManager("count")-1);
		roiManager("rename", "Cell"+num_cell+1+"_spot"+i_r+1+"_trans");
	}

	if( th_blk != -1 ) {
		int_spot = newArray(nResults); 
		selectWindow("MinProj");
		for(i_r = nResults-1; i_r >=0; i_r--){
			if( getPixel(getResult("X",i_r),getResult("Y",i_r)) > th_blk ){
				roiManager("select", i_r+1); // gap between spots and results: cell ROI
				roiManager("delete");
			}
		}
	}

	// Offers to the user the possibility to remove spots
	selectWindow(tit_img);
	waitForUser("Remove spots you do not want to analyze");
	run("Clear Results");
	roiManager("deselect");
	roiManager("measure");

	// updates results table and roiManager according to user choice
	coordinates = newArray(3*(nResults-1));
	for(i_r = 0; i_r < nResults-1; i_r++){ // first ROI is the cell
		coordinates[i_r*3] = getResult("X",i_r+1);
		coordinates[i_r*3+1] = getResult("Y",i_r+1);
		selectWindow("Assoc_slices_trans");
		val_z = getPixel(coordinates[i_r*3],coordinates[i_r*3+1]);
		coordinates[i_r*3+2] = val_z;
		roiManager("select", i_r+1);
		roiManager("rename", "Cell"+num_cell+1+"_spot"+i_r+1+"_trans");
	}	
	return coordinates; // (X,Y,Z) coordinates, every 3 indexes is a new point
}

// function to detect the pigments based on FindMaxima function on the Max projection
// Spots with intensity bellow the mean value of all spots are removed
function detectPigmentRedCh(tit_img,tit_img_Zinfo,prom_param,prox_th_XY,width, height){

	selectWindow("MaxProj");
	roiManager("select", 0);
	run("Find Maxima...", "prominence="+prom_param+" output=List");

	// search for the z-value corresponding to detected spots
	selectWindow(tit_img_Zinfo);
	for(i_r = 0; i_r < nResults; i_r++){
		val_z = getPixel(getResult("X",i_r),getResult("Y",i_r));
		setResult("Z", i_r, val_z);
	}
	updateResults();
}

// find for each pixel in tit_img_min which slice of tit_img contains the value
// ---> gives the z-slice of the min/max projection
function findCorrespondanceProjImage(tit_img,tit_img_proj,width, height,slices,name_assoc_img){
	tab_min = newArray(width*height);
	tab_index = newArray(width*height);
	// reads the values of the projection image
	selectWindow(tit_img_proj);
	for(i_w = 0; i_w < width; i_w++){
		for(i_h = 0; i_h < height; i_h++){
			tab_min[i_w*height+i_h] = getPixel(i_w,i_h);
		}
	}
	// looks the value for each pixel on each slice and compare it to the projection
	selectWindow(tit_img);
	for(i_s=1; i_s <= slices; i_s++){
		Stack.setSlice(i_s);
		for(i_w = 0; i_w < width; i_w++){
			for(i_h = 0; i_h < height; i_h++){
				if(getPixel(i_w,i_h) == tab_min[i_w*height+i_h])
					tab_index[i_w*height+i_h] = i_s; 
			}
		}
	}
	// creates the image with the associated slice number for each pixel of the projection
	newImage(name_assoc_img, "8-bit black", width, height, 1);
	for(i_w = 0; i_w < width; i_w++){
		for(i_h = 0; i_h < height; i_h++){
			setPixel(i_w, i_h, tab_index[i_w*height+i_h]);
		}
	}
}

// function to display the results for the concerned color: max/mean/positive/proportion of positive pixels
// we consider the pigment is positive for this channel if its maximum is above the threhsold stored in variable
// max_bg_color
function displayResults(str_color,res_line,max_value_color,mean_value_color,max_bg_color,prop_pos_pix_color){
	setResult("Max "+str_color+" around spot",res_line,max_value_color);
	setResult("Mean "+str_color+" around spot",res_line,mean_value_color);
	if( max_value_color>max_bg_color )
		setResult(str_color+" positive",res_line,"yes");
	else 
		setResult(str_color+" positive",res_line,"no");
	setResult("Proportion "+str_color+" positive",res_line,prop_pos_pix_color);
}

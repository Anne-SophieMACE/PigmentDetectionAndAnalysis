/*
 * macro distancePigments_nucleus_v1.ijm
 * 
 * This macro computes the distance between pigments and nucleus
 * pigments already detected thanks to macro detectPigment_quantifySignalAround_vf.ijm
 * -> ROIs (cells and pigments) are stored in a zip file that should be
 * loaded by the user at the beginning of the macro
 * 
 * The last channel should be the nucleus one (= the cyan)
 *
 */

// clears everything
run("Close All");
run("Clear Results"); 

run("Set Measurements...", "centroid stack redirect=None decimal=2");
// asks the user the ROI file he wants to analyze
file_ROI = File.openDialog("Choose your file of spots ROIs"); 
if( endsWith(file_ROI, ".zip") == false )
	exit("You should choose a .zip file containing ROIs");
	
// resets manager
roiManager("reset");
roiManager("Show All");
roiManager("Show None");

// Dialog box for the parameters (steps in xy and z)
Dialog.create("Voxel size");
Dialog.addNumber("Pixel size (µm):", 0.0645);
Dialog.addNumber("Z-step (µm):", 0.2);
Dialog.show();
xy_step = Dialog.getNumber();
z_step = Dialog.getNumber();

// opens the ROI file
roiManager("open", file_ROI);

// opens the image file corresponing to the acquisiton
orig_img = substring(file_ROI,lastIndexOf(file_ROI, File.separator)+1,indexOf(file_ROI, "_ROI_"));
open(orig_img);
orig_img_name = getTitle();
img_path = getDirectory("image");
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // scale in pixels
getDimensions(width, height, channels, slices, frames); // image dimensions

// check if there is a "slice" information 
str_slice = "1-"+slices;
if( indexOf(file_ROI, "_slices") != -1 ){ // slice information -> substack made
	str_slice = substring(file_ROI,indexOf(file_ROI, "_slices")+lengthOf("_slices"),lastIndexOf(file_ROI, "."));
}
// duplicates the slices
selectWindow(orig_img_name);
run("Select All");
run("Duplicate...", "title=img_analysis duplicate slices="+str_slice);
img_name = getTitle();

// deletes the original image
selectWindow(orig_img_name);
close();

nbROIbefore_draw = roiManager("count");

// searches for the number of cells i.e. the number of ROIs 
// with name starting by "ROI"
nbROI_cell = 1;
if( indexOf(file_ROI, "cell") != -1 ){
	roiManager("select", nbROI_cell);
	roi_name = Roi.getName();
	while( startsWith(roi_name,"ROI") ){
		nbROI_cell++;
		roiManager("select", nbROI_cell);
		roi_name = Roi.getName();
	}
}

// removing the first ROIs (= ROIs of the cells of the first macro)
for( i_r = 0; i_r < nbROI_cell; i_r++){ 
	roiManager("select", 0);
	roiManager("delete");
}

// measurement of the spot = pigments coordinates
roiManager("deselect");
roiManager("measure");

nbPigments = roiManager("count");
spots_xcoord = newArray(nbPigments);
spots_ycoord = newArray(nbPigments);
spots_zcoord = newArray(nbPigments);
spot_roi_name = newArray(nbPigments); // to know which cell it corresponds to

// get coordinates of the pigment centers
for (i_r = 0; i_r < nbPigments; i_r++) {
	// X,Y are given in pixels (because unscaled image)
	spots_xcoord[i_r] = getResult("X", i_r)*xy_step;
	spots_ycoord[i_r] = getResult("Y", i_r)*xy_step;
	roiManager("select", i_r);
	spot_roi_name[i_r] = Roi.getName;
	// Slice is given in numero of slice
	spots_zcoord[i_r] = getResult("Slice", i_r)*z_step;
}

// re-opens ROI zip to keep only cells ROIs
roiManager("reset");
roiManager("open", file_ROI);
while( roiManager("count")!= nbROI_cell){
	roiManager("select", roiManager("count")-1);
	roiManager("delete");
}

// detect Nuclei in the DAPI
selectWindow(img_name);
roiManager("show all");
roiManager("show none");
run("Duplicate...", "title=nucl_img duplicate channels="+channels);
// creation of the mask image
run("Z Project...", "projection=[Max Intensity]");
rename("mask_nucl");
run("Subtract Background...", "rolling=500");
setOption("ScaleConversions", true);
run("8-bit");
run("Auto Threshold", "method=Otsu white");
run("Remove Outliers...", "radius=2 threshold=50 which=Bright");

// table of coordinates for the nuclei
nucleus_coord_X = newArray(nbROI_cell);
nucleus_coord_Y = newArray(nbROI_cell);
nucleus_coord_Z = newArray(nbROI_cell);
nucleus_majorAxis = newArray(nbROI_cell);
nucleus_minorAxis = newArray(nbROI_cell);
for (i_roi = 0; i_roi < nbROI_cell; i_roi++) {
	run("Clear Results"); 
	// look for one nucleus in the cell ROI
	selectWindow("mask_nucl");
	roiManager("select", i_roi);
	run("Analyze Particles...", "size=100-Infinity add");
	if( roiManager("count") != nbROI_cell+1 ){ // several or none nucleus found
		while( roiManager("count")!= nbROI_cell ){ // delete all ROIs found by Analyze particles
			roiManager("select", roiManager("count")-1);
			roiManager("delete");
		}
		selectWindow("nucl_img");
		waitForUser("Draw your nucleus (at the correct Z), add it in the manager");
		while( roiManager("count") != nbROI_cell+1 ){
			waitForUser("Draw your nucleus (Cell"+i_roi+1+") (at the correct Z), add it in the manager");
		}
	}
	else{ // only one ROI, we suppose it is correct
		// the slice with bigger std is supposed to be the focal plan
		run("Set Measurements...", "standard redirect=None decimal=2");
		selectWindow("nucl_img");
		run("Measure Stack...");
		std_tab_max = getResult("StdDev", 0);
		slice_std_max = getResult("Slice", 0);
		for (i_res = 1; i_res < nResults; i_res++) {
			if( getResult("StdDev", i_res) > std_tab_max ){
				std_tab_max = getResult("StdDev", i_res);
				slice_std_max = getResult("Slice", i_res);
			}
		}
		selectWindow("nucl_img");
		Stack.setSlice(1);
		roiManager("show all");
		// change the name so that it integrates the slice
		roiManager("select", roiManager("count")-1);
		if( lengthOf(toString(slice_std_max)) == 1 )
			roiManager("rename", "000"+slice_std_max+"-"+Roi.getName);
		else
			roiManager("rename", "00"+slice_std_max+"-"+Roi.getName);
		// asks the user to check
		roiManager("select",roiManager("count")-1);
		waitForUser("Check the nucleus (Cell"+i_roi+1+"): if uncorrect draw the correct one (at the correct Z), add it in the manager");
	}

	run("Clear Results"); 
	roiManager("select", roiManager("count")-1); // nucleus
	run("Set Measurements...", "centroid fit stack redirect=None decimal=2");
	roiManager("measure");

	// reads all results and converts them to µm
	nucleus_coord_X[i_roi] = getResult("X", 0)*xy_step;
	nucleus_coord_Y[i_roi] = getResult("Y", 0)*xy_step;
	nucleus_coord_Z[i_roi] = getResult("Slice", 0)*z_step;
	nucleus_majorAxis[i_roi] = getResult("Major", 0)*xy_step;
	nucleus_minorAxis[i_roi] = getResult("Minor", 0)*xy_step;
	
	roiManager("deselect");
	// removes the nucleus but saves the nuclei ROI
	roiManager("select", roiManager("count")-1);
	roiManager("save", img_path+File.separator+orig_img+"_nucleus"+i_roi+1+".roi");
	while( roiManager("count")!= nbROI_cell ){
		roiManager("select", roiManager("count")-1);
		roiManager("delete");
	}
}

// fills the result table of distances
run("Clear Results");
for (i_r = 0; i_r < nbPigments; i_r++) {
	cell_roi_number = parseInt(substring(spot_roi_name[i_r], 4, indexOf(spot_roi_name[i_r], "_")))-1;
	setResult("Roi name", i_r, spot_roi_name[i_r]);
	setResult("Distance to nucleus (µm)", i_r, comp_dist(nucleus_coord_X[cell_roi_number],nucleus_coord_Y[cell_roi_number],nucleus_coord_Z[cell_roi_number],spots_xcoord[i_r],spots_ycoord[i_r],spots_zcoord[i_r]));
}
// save pigments results
saveAs("Results", img_path+File.separator+orig_img+"distanceNucleus_eachPigment.xls");

// fills the result table for nucleus information
run("Clear Results");
for (i_nucl = 0; i_nucl < nbROI_cell; i_nucl++) {
	setResult("X coord (µm)", i_nucl, nucleus_coord_X[i_nucl]);
	setResult("Y coord (µm)", i_nucl, nucleus_coord_Y[i_nucl]);
	setResult("Z coord (µm)", i_nucl, nucleus_coord_Z[i_nucl]); 
	setResult("Major Axis (µm)", i_nucl, nucleus_majorAxis[i_nucl]);
	setResult("MinorAxis (µm)", i_nucl, nucleus_minorAxis[i_nucl]);
	setResult("Approx. radius (µm)", i_nucl, (nucleus_minorAxis[i_nucl]+nucleus_majorAxis[i_nucl])/4);

	roiManager("open", img_path+File.separator+orig_img+"_nucleus"+i_nucl+1+".roi");
	roiManager("select", roiManager("count")-1);
	roiManager("rename", "Nucleus"+i_nucl+1);
	del = File.delete(img_path+File.separator+orig_img+"_nucleus"+i_nucl+1+".roi");
}
// saves ROIs and results
roiManager("save", img_path+File.separator+orig_img+"_cells_nuclei.zip")
saveAs("Results", img_path+File.separator+orig_img+"nucleus_information.xls");

// function that computes distance between two points of coordinates 
// (x1,y1,z1) and (x2,y2,z2)
function comp_dist(x1,y1,z1,x2,y2,z2){
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}

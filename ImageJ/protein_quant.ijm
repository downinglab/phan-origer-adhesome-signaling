// --- ImageJ/FIJI Macro: Consolidated 4-Channel Analysis with Group Tracking ---

// Setup Input and Output Directories
inputDir = getDirectory("Choose your INPUT folder containing .oir files");
outputDir = getDirectory("Choose your OUTPUT folder for results");
list = getFileList(inputDir);

// Configure global measurement settings
run("Set Measurements...", "area mean integrated redirect=None decimal=3");
run("Clear Results");

setBatchMode(true); // Hides windows during processing

// Loop through all files in the folder
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".oir")) {
        
        showProgress(i+1, list.length);
        
        // Open file using Bio-Formats
        path = inputDir + list[i];
        run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
        
        baseName = getTitle();
        
        // Figure out condition from filename
        condition = "Unknown";
        if (indexOf(baseName, "LZ") >= 0) condition = "lacZ KD";
        if (indexOf(baseName, "S3") >= 0) condition = "SHROOM3 KD";
        
        // Dynamic Group Extraction
        groupLabel = "Unknown";
        gIndex = indexOf(baseName, "-g");
        if (gIndex >= 0) {
            subStr = substring(baseName, gIndex + 1); 
            nextDash = indexOf(subStr, "-");
            if (nextDash >= 0) {
                groupLabel = substring(subStr, 0, nextDash);
            }
        }
        
        // Split Channels
        run("Split Channels");
        
        // Assign windows based on original 4-channel metadata mapping
        wDAPI  = "C1-" + baseName;
        wDesmo = "C2-" + baseName;
        wSem7A = "C3-" + baseName;
        wFibro = "C4-" + baseName;
        
        // --- Step 4: Create the Expanded DAPI Mask ---
        selectWindow(wDAPI);
        run("Duplicate...", "title=DAPI_Mask");
        run("Gaussian Blur...", "sigma=2");
        
        // Automatically threshold the DAPI
        setAutoThreshold("Otsu dark");
        run("Convert to Mask");
        
        // Expand the mask into the cytoplasm (Dilate 4 times)
        for (d = 0; d < 4; d++) {
            run("Dilate");
        }
        
        // Create selection from the expanded mask and add to ROI Manager
        run("Create Selection");
        roiManager("add");
        
        // --- Step 5: Measure Protein Channels ---
        channels = newArray(wEphA2, wEphA7, wFN1);
        proteins = newArray("EphA2", "EphA7", "FN1");
        
        for (c = 0; c < channels.length; c++) {
            selectWindow(channels[c]);
            roiManager("select", 0);
            run("Measure"); // Adds data to the native Results table safely
            
            // Appends custom tracking variables to the active row
            currentRow = nResults - 1;
            setResult("image", currentRow, list[i]);
            setResult("condition", currentRow, condition);
            setResult("group", currentRow, groupLabel);
            setResult("target", currentRow, proteins[c]);
        }
        
        // Clean up windows and ROI manager for the next image loop
        roiManager("reset");
        run("Close All");
    }
}

// Save the standard results table cleanly as your CSV
saveAs("Results", outputDir + "cell-signaling-if.csv");

setBatchMode(false);
print("Finished! All images combined into: cell-signaling-if.csv");
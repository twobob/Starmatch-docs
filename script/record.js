var SAVED_NUMBER_PRECISION = 5; // old value 5
var STORED_VALUES_INTERNAL_PRECISION = 5; //old value 5
var TARGET_PROCESSED_DATA_PRECISION = 4; // old value 4
var POSITIONS_PRECISION = 5; //old value 5
var THEMES_PRECISION = 4; //old value 4
var FUNCTION_PRECISION = 5; // d.p. in fn coefficient values
var OPTIONS_ARRAY_SIZE = 13; // was 13
var DEBUG = 1; // 0 shuts it up
var curveDataSet = []; // GLOBAL STORE FOR COEFFICIENTS
var inflectionPoints = [];
var SubjectValues = [];
var NonRadioOptionsObject = {
}; // Defaults for orbType 1 (planet orbs)
var OptionsObject = {
    aoIndex: 0,
    tfIndex: 0,
    precessionFlag: 0,
    orbType: 0,
    poIndex: 0,
    chartType: 0,
}; // Defaults for orbType 1 (planet orbs)
/* poIndex = 0, 1 (Lilly, al-Biruni) */
var OptionObjectValues = {
    aoIndex: [0, 1, 2, 3, 4],
    tfIndex: [0, 4, 8],
    precessionFlag: [0, 1],
    orbType: [0, 1],
    poIndex: [0, 1],
    chartType: [0, 1],
};
var OptionsObjectLon = 0;
var OptionsObjectLat = 0;
var OptionsObjectDateTime = 0;
var OptionsObjectDateTimeAccurate = 0;
var OptionsObjectUtOffset = 0;
var OptionsObjectDstOffset = 0;
var GlobalCurrentSaveOptions = [];
var delimeter = "¬"; // IS GLOBAL 
const __PLAYLIST__PREFIX__ = "_æßÇ_"; // converted to escape chars below. However I need to test this before making any adjustments.
const __PREFIX__ = "\u00bb"; // <-- STARMATCH MOBILE NAMESPACE.
var UniqueSubjectProcessedData = [];
var selectionsList = [];
var lastFullChartReport = [];
function _f__TheCurrentResearchOptionsMatchedTheStoredOnes(optionsArray) {
    var offset = 6; // What's this?
    var AllTheOptionsMatch = true;
    var amountOfResearchOptions = Object.size(optionsArray) - offset;
    for (z = 0; z < amountOfResearchOptions; z++) {
        if (parseFloat(GlobalCurrentSaveOptions[z + offset]) != parseFloat(optionsArray[z + offset])) {
            AllTheOptionsMatch = false;
        }
    }
    return AllTheOptionsMatch;
}
var currentSelectionArray = [];

function LoadChosenSelection(nameOfSelectionList) {
    if (typeof nameOfSelectionList === "undefined") {
        _f__infoP("enter name of selection group in the text box");
        return;
    }
    $('#selection_name').val(nameOfSelectionList);
    _f__populateSelectionsList();
    $('#filterTxt').val('');
    _f__populateRecordsList();
    var nameOfSelectionGroup = __PLAYLIST__PREFIX__ + nameOfSelectionList;
    var nameOfRecord = nameOfSelectionList;
    if (localStorage.getItem(nameOfSelectionGroup) !== null) {
        var fields = localStorage.getItem(nameOfSelectionGroup).split(delimeter);
        $('#person option').prop('selected', false);
        jQuery.each(fields, function (i, field) {
            $('#person option').filter(function (index, e) {
                return $(e).text() == decodeURIComponent(field)
            }).prop('selected', true);
        });
    }
}

function ClearCurrentSelection() {
    $(".personRecords option:selected").prop("selected", false);
    $("#selections li").removeClass("selected");
}

function SaveCurrentSelection() {
    if ($('#selection_name').val() === "undefined") {
        _f__warnP("Enter a name for the selection group");
        return;
    }
    if ($('#selection_name').val().trim() == "") {
        _f__warnP("Enter a name for the selection group");
        return;
    }
    var count = $('#person  :selected').length
    var length = $('#person > option').length;
    if (count < 1) {
        _f__warnP("Select some records to save as a selection group");
        return;
    }
    var currentSelectionArray = $('#person  :selected');
    var nameOfSelectionGroup = __PLAYLIST__PREFIX__ + $('#selection_name').val();
    var nameOfRecord = $('#selection_name').val();
    if (localStorage.getItem(nameOfSelectionGroup) !== null) {
        _f__infoP("Updated selection group: " + nameOfRecord, 3000);
        _f__deleteRecordFromSelectionList(nameOfSelectionGroup);
    } else {
        _f__infoP("Saving selection group: " + nameOfRecord, 3000);
    }
    var exportString = "";
    $.each(currentSelectionArray, function (index, thing) {
            exportString += encodeURIComponent(thing.innerHTML) + delimeter;
        })
    if (exportString.length > 1) {
        exportString = exportString.slice(0, -1);
    } else {
    }
    currentSelectionArray = exportString;
    localStorage.setItem(nameOfSelectionGroup, currentSelectionArray);
    _f__populateSelectionsList();
}

function _f__SelectAllRecords() {
    $("#selectionsHolder").dialog("open");
}

function DeleteCurrentSelection() {
        if (!$('#selections li').hasClass("selected")) {
            _f__warnP('Select a record to delete ');
        }
        var thingToDel = $('#selections li.selected').text();
        if (localStorage.getItem(__PLAYLIST__PREFIX__ + thingToDel)) {
            localStorage.removeItem(__PLAYLIST__PREFIX__ + thingToDel)
        };
        _f__populateSelectionsList();
    }
function _f__populateSelectionsList() {
    selectionsList = document.getElementById('selections');
    var usePreviousSelection = false;
    var lastSelectedVal = "";
    if ($('#selections li').hasClass("selected")) {
        usePreviousSelection = true;
        lastSelectedVal = $('#selections li.selected').text();
    }
    $('#selections').empty();
    var name = new Array();
    if (_f__storageAvailable('localStorage') == false) {
        _f__warnP('Too bad, no localStorage for us. ');
        return;
    }
    for (var i = 0; i <= localStorage.length - 1; i++) {
        var key = localStorage.key(i);
        var val = localStorage.getItem(key);
        if (key.includes(__PLAYLIST__PREFIX__)) {
            var tidyName = key.replace(__PLAYLIST__PREFIX__, '');
            var opt = document.createElement('li');
            if (lastSelectedVal == tidyName) {
                opt.className = "selected inlinePersonRecordClass";
            } else {
                opt.className = "inlinePersonRecordClass";
            }
            opt.value = tidyName;
            opt.onclick = function () {
                $('#selections li').removeClass('selected');
                $(this).addClass('selected');
                LoadChosenSelection($(this).text())
            };
            opt.title = "Select the group: " + tidyName;
            opt.innerHTML = tidyName;
            selectionsList.appendChild(opt);
        }
    }
}

function _f__updateThemeValuesSHUNT() {
    _f__infoP("Mass Update In Progress", 2000);
    var timeoutID = window.setTimeout(_f__updateThemeValues, 300);
}

function _f__updateThemeValues() {
    _f__SetOptionsObjectFromWidgetState();
    var subjectFlag = 0; // global hack to indicate 1st record in list (the subject)
    recordsList = document.getElementById('person');
    var index = $('#person option'); // all of them
    var currentSubjectName = $('nav select option').first().text().toString();
    index.each(function (innerIndex, value) {
        var $this = $(this);
        currentRecordName = value.text;
        var recordHolder = _f__decode(localStorage[currentRecordName]).split(delimeter).slice(24)
            /*  recordHolder is a slice of the current records current LON LAT D/TOB TIMEACC UTOFFSET DSTOFFSET and the "potentially any length" options object flags
0: "-1.132"  LON
1: "53.522"  LAT
2: "0" D/TOB
3: "0" ACC
4: "0" UTOFFSET
5: "0" DSTOFFSET
6: "0" 
7: "0"
8: "0"
9: "0"
10: "0"
11: "0"
12: "0"
length: 13	// now OPTIONS_ARRAY_SIZE

            var passedtimeAccurate = optionsArray[3];
            var passedlatitude = optionsArray[1];
            var passedlongitude = optionsArray[0];

            var passedUt = optionsArray[4];
            var passedDst = optionsArray[5];

            var passedaoIndex = optionsArray[6];
            var passedtfIndex = optionsArray[7];
            var passedprecessionFlag = optionsArray[8];

var passedorbType = optionsArray[9];	


	*/
        var profileValues = [];
        var positions = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        var LocalTheme = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        var optionsArray = new Array(OPTIONS_ARRAY_SIZE); // num options?
        var CurrentTimeOfBirthAccurate = false;
        _f__StorageRetrieveOperation(value.text, positions, LocalTheme, optionsArray);
        CurrentTimeOfBirthAccurate = _f__BirthTimeAccurate(value.text);
        if (!_f__RecordEntriesAreNoneZero(positions)) {
			
			
			
            _f__infoP("I Booted the record " + value.text + " It has no meaningful POSITIONS data ", 10000);
        }
        else {
            optionsArray[6] = OptionsObject.aoIndex;
            optionsArray[7] = OptionsObject.tfIndex;
            optionsArray[8] = OptionsObject.precessionFlag;
            optionsArray[9] = OptionsObject.orbType;
            optionsArray[10] = OptionsObject.poIndex;
            optionsArray[11] = OptionsObject.chartType;
            _f__SetNativityYearAndOptionsObject(value.text);
            getThemeValues(positions[0], positions[1], positions[2], positions[3], positions[4], positions[5], positions[6], positions[7], positions[8], positions[9], positions[10], positions[11]);
            localStorage.removeItem(currentRecordName)
            _f__StorageSaveOperation(currentRecordName, positions, theme, optionsArray);
        }
    });
    _f__infoP("MASS UPDATE complete", 1500);
}

function SaveCurrentSelectionWithName(NameToSave) {
    var currentSelectionArray = $('#person  :selected');
    var nameOfSelectionGroup = __PLAYLIST__PREFIX__ + NameToSave;
    var nameOfRecord = NameToSave;
    if (localStorage.getItem(nameOfSelectionGroup) !== null) {
        _f__infoP("Updated selection group: " + nameOfRecord, 3000);
        _f__deleteRecordFromSelectionList(nameOfSelectionGroup);
    } else {
        _f__infoP("Saving selection group: " + nameOfRecord, 3000);
    }
    var exportString = "";
    $.each(currentSelectionArray, function (index, thing) {
            exportString += encodeURIComponent(thing.innerHTML) + delimeter;
        })
    if (exportString.length > 1) {
        exportString = exportString.slice(0, -1);
    } else {
    }
    currentSelectionArray = exportString;
    localStorage.setItem(nameOfSelectionGroup, currentSelectionArray);
    _f__populateSelectionsList();
}

function timeStamp() {
    var now = new Date();
    var date = [now.getMonth() + 1, now.getDate(), now.getFullYear()];
    var time = [now.getHours(), now.getMinutes(), now.getSeconds()];
    var suffix = (time[0] < 12) ? "AM" : "PM";
    time[0] = (time[0] < 12) ? time[0] : time[0] - 12;
    time[0] = time[0] || 12;
    for (var i = 1; i < 3; i++) {
        if (time[i] < 10) {
            time[i] = "0" + time[i];
        }
    }
    return date.join("_") + "_" + time.join("-") + "_" + suffix + ".txt";
}
var TotalIsInLimits = 0;
var TotalTotalTests = 0;
var TotalOutputString = "";
var Records = [];

function _f__RunAllTestsSHUNT() {
    _f__infoP("Beginning Tests", 2000);
    var timeoutID = window.setTimeout(_f__RunAllTests, 300);
}

function _f__RunAllTests() {
	
		if (OptionsObject.chartType != 0)
		{
				_f__warnP("SET THE CHART TYPE TO 0", 10000);
				return;
		}

	resizeDoneOnce = false;
	
        TotalIsInLimits = 0;
        TotalTotalTests = 0;
        TotalOutputString = "";
        _f__clearCharts();
        _f__debugClear(); // Do NOT clear debug proactively. WJ18 request 31/01/06 (reinstated 22/03/2018 - Will)
        _f__warnClear();
        var RecordToSort = {
            'sortValuesWithString': [],
            'header': ""
        };
        Records = [];
        var holder = [];
		
        var parsedName = "";
        $('.inlinePersonRecordClass').each(
            function (inx, thing) {
                parsedName = thing.title.replace('Select the group: ', '');
                if (parsedName.endsWith('TEST')  && ($('#runTEST').is(':checked'))    ) {
                    holder[holder.length] = [parsedName.replace('_F_R','').replace('_M_R','').replace('_F','').replace('_M','').replace('TEST', '').trim(), parsedName, "ALL"]
					holder[holder.length-1].SelectionGroup = parsedName.replace('TEST').trim();
                };
                if (parsedName.endsWith('TESTM')  && ($('#runTESTM').is(':checked'))) {
                    holder[holder.length] = [parsedName.replace('_F_R','').replace('_M_R','').replace('_F','').replace('_M','').replace('TESTM', '').trim(), parsedName, "ALL_F"]
					holder[holder.length-1].SelectionGroup = parsedName.replace('TESTM').trim();
                };
                if (parsedName.endsWith('TESTF') && ($('#runTESTF').is(':checked'))) {
                    holder[holder.length] = [parsedName.replace('_F_R','').replace('_M_R','').replace('_F','').replace('_M','').replace('TESTF', '').trim(), parsedName, "ALL_M"]
					holder[holder.length-1].SelectionGroup = parsedName.replace('TESTF').trim();
                };
            });
			
			debugP("Beginning "+holder.length+ " tests");
        for (var l = 0; l < holder.length; l++) {
			
			
			if (localStorage.getItem(holder[l][0]) == null)
			{
				_f__warnP("Test for "+holder[l][0]+" Has no stored record with that name", 10000);
			}
			 var CheckTimeOfBirthAccurate = true;
            			
				
                 if ($('#TimeImportant').prop("checked")){
	CheckTimeOfBirthAccurate = _f__BirthTimeAccurate(holder[l][0]);


			 if (CheckTimeOfBirthAccurate  != 1){
			
continue;


			
	}
	
/*
				 if( CheckTimeOfBirthAccurate == false) { 
				 debugP("Skipping record " + holder[l][0] + " for inaccurate TOB");
				  debugP("------------------------");
                    continue;
				 }
*/
				 }

            $('#SubjectName').val(holder[l][0]);
            $('#selections li').removeClass('selected');
            $('#person option').prop('selected', false);
            $('#person option').filter(function (index, e) {
                return $(e).text() == holder[l][0]
            }).prop('selected', true);
            _f__ClickedLoadForTests('subject');
            LoadChosenSelection(holder[l][1]);
			
			
		    
			_f__ClickedxProfileMean(holder[l][2]);
        
		if (TotalIsInLimits > 0){
			debugP (holder[l].SelectionGroup.replace('undefined','')  );
        }
		
		}
        TotalOutputString = "";
        TotalOutputString += TotalIsInLimits + " Matches out of " + TotalTotalTests + " tests total " + '\n';
        debugP(TotalIsInLimits + " Matches out of " + TotalTotalTests + " tests total ");
	
	for (t = 0; t < Records.length; t++) {

			if (Records[t].sortValuesWithString.length<1)
			continue;
	
			
			TotalOutputString += Records[t].header;
            
            var myOutput = Records[t].sortValuesWithString.sort(function (a, b) {
                return a[0] < b[0];
            });
           
            for (a = 0; a < myOutput.length; a++) {
				
                TotalOutputString += myOutput[a][1];
            
			}
        }
        ChosenFilename = "ALL_TESTS.txt";
		
		
       
        _f__infoP(TotalIsInLimits + " Matches out of " + TotalTotalTests + " tests total ", 5000);
		 _f__createDownloadableFile(ChosenFilename, TotalOutputString);
        _f__ResetTheEngine();
		
    }
   
function _f__ClickedxProfileMean(ListToCall) {
debugP("I CALLED _f__ClickedxProfileMean method"+method);
		
		
		var reservedPrintout=[];
		
        var RecordToSort = {
            'sortValuesWithString': [],
            'header': ""
        };
        var StoredSelection = $("#selection_name").val();
        if (!($('nav select option:selected').length))
            return;
       
        var ChosenFilename = "PT_" + $('#SubjectName').val() + "_" + $("#selection_name").val() + "_" + timeStamp();
       
        if (!($('nav select option:selected').length))
            return;
        var subjectFlag = 0; // global hack to indicate 1st record in list (the subject)
        if (chartType == 1) // display purposes only, type 0 needed here
        {
            _f__warnP("select chart type 0 in Options (ALT+SHIFT+o)");
         
            return;
        }
        if (typeof (SubjectProcessedData) === "undefined" || !_f__RecordEntriesAreNoneZero(SubjectsValues)) {
            _f__warnP("Load a SUBJECT to profile");
            return;
        }
        recordsList = document.getElementById('person');
        var index = $('#person option:selected');
        if (index.length < 1) {
            _f__warnClear();
            _f__warnP("Please select TARGET record(s) from the list to cross-profile");
            return;
        }
        _f__CreateOptionsObject();
        var fitArray = [];
        var profileValues = [];
        _f__DoCompleteOptionsSaveThenLoadCycle();
        var selectedIndexLength = index.length;
        var currentSubjectName = $('#SubjectName').val().toString();
        var currentRecordName = "";
        var reprocessedRecords = [];
        var totalFail = 0;
        var br = "<br />";
        var subReporter = "Chart for " + currentSubjectName + CRLF;
        var reporter = [];
        reporter[reporter.length] = subReporter;
        subReporter = "";
        index.each(function (innerIndex, value) {
            var $this = $(this);
            currentRecordName = value.text;
            if (currentSubjectName == currentRecordName) {
                if (selectedIndexLength < 2) {
                    totalFail += 1;
                }
                return;
            }
            _f__LoadSubjectChartData(currentSubjectName);
            if ($this.length) {
                var positions = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                var LocalTheme = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                var optionsArray = new Array(OPTIONS_ARRAY_SIZE); // num options?
                var CurrentTimeOfBirthAccurate = false;
                _f__StorageRetrieveOperation(value.text, positions, LocalTheme, optionsArray);
                CurrentTimeOfBirthAccurate = _f__BirthTimeAccurate(value.text);
                if ($('#TimeImportant').prop("checked") && CurrentTimeOfBirthAccurate == false) { //debugP("Skipping record " + value.text + " for inaccurate TOB");
                    return;
                }
                if (!_f__RecordEntriesAreNoneZero(positions)) {
                    _f__infoP("I Booted the record " + value.text + " It has no meaningful POSITIONS data ", 10000);
                }
                else {
                    _f__SetNativityYearAndOptionsObject(value.text);
                    if (!_f__TheCurrentResearchOptionsMatchedTheStoredOnes(optionsArray)) {
                        /*  _f__infoP("Used generated theme data for "+currentRecordName); */
                        /* 	//debugP("Used generated Path!"); */
                        getThemeValues(positions[0], positions[1], positions[2], positions[3], positions[4], positions[5], positions[6], positions[7], positions[8], positions[9], positions[10], positions[11]);
                        LocalTheme = theme.slice();
                        _f__CreateOptionsObjectWithPerRecordDataIntact(optionsArray);
                    }
                    var sArrayTcopy = [
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], SubjectProcessedData
                    ];
                    var tArrayTcopy = [
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], LocalTheme
                    ];
debugP("record.js are we here? 1");	// not on xprofile button call...
                    if (method < 3) // 2 and 3 are external to engine and control graph display
                    {
                        var profileData = xProfile(sArrayTcopy, tArrayTcopy);
                        var barColour = profileData.pCoincidence;
                        if (!subjectFlag) {
                            inflectionPoints.push(profileData.sInflections)
                            inflectionPoints.push(profileData.sArrayTp);
                            inflectionPoints.push(profileData.sArray);
                        }
                        inflectionPoints.push(profileData.tInflections)
                        inflectionPoints.push(profileData.tArrayTp);
                        inflectionPoints.push(profileData.tArray);
                        subjectFlag = 1; // hack to ensure subject only pushed once				
                        var numberToSave = profileData.fitC + profileData.fitS; // this is the actual 'fit'
                        if ((profileData.fitC == 1) && (profileData.fitS == -1)) {
                            debugP("Error: " + currentRecordName + " no corresponding theme values outside +/- std. dev.");
                            debugP("See type 1 charts: unable to calculate fit value");
                        } else { // do not put error case names on list, pass btAccurate for graph colour
                            profileValues.push([currentRecordName, numberToSave, barColour, CurrentTimeOfBirthAccurate.toString()]);
                        }
                        subReporter = _f__pad(currentRecordName, 30, " ") + _f__pad(numberToSave.toPrecision(SAVED_NUMBER_PRECISION), 10, " ");
                        reporter[reporter.length] = subReporter;
                    }
					
					
                }
            }
        });
        if (totalFail > 0) {
            return;
        }
        profileValues.sort(function (a, b) {
            return b[1] - a[1];
        });
        for (n = 0; n < profileValues.length; n++)
            fitArray.push(profileValues[n][1]);
        if (method < 3) { // allow meaninglessly small data sets e.g. 1 item, no warning
            var storeFirstData = [];
            for (var i = 0; i < profileValues.length; i++) {
                storeFirstData[storeFirstData.length] = [profileValues[i][0], profileValues[i][1]];
            }
        }
        LoadChosenSelection(ListToCall);
        var subjectFlag = 0; // global hack to indicate 1st record in list (the subject)
        recordsList = document.getElementById('person');
        var index = $('#person option:selected');
        _f__CreateOptionsObject();
        var fitArray = [];
        var profileValues = [];
       
        _f__DoCompleteOptionsSaveThenLoadCycle();
        var selectedIndexLength = index.length;
        var currentSubjectName = $('#SubjectName').val().toString();
        var currentRecordName = "";
        var reprocessedRecords = [];
        var totalFail = 0;
        var br = "<br />";
        var subReporter = "Chart for " + currentSubjectName + CRLF;
        var reporter = [];
        reporter[reporter.length] = subReporter;
        subReporter = "";
        index.each(function (innerIndex, value) {
            var $this = $(this);
            currentRecordName = value.text;
            _f__LoadSubjectChartData(currentSubjectName);
            if ($this.length) {
                var positions = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                var LocalTheme = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                var optionsArray = new Array(OPTIONS_ARRAY_SIZE); // num options?
                var CurrentTimeOfBirthAccurate = false;
                _f__StorageRetrieveOperation(value.text, positions, LocalTheme, optionsArray);
                CurrentTimeOfBirthAccurate = _f__BirthTimeAccurate(value.text);
                if ($('#TimeImportant').prop("checked") && CurrentTimeOfBirthAccurate == false) {
                    return;
                }
                if (!_f__RecordEntriesAreNoneZero(positions)) {
                    _f__infoP("I Booted the record " + value.text + " It has no meaningful POSITIONS data ", 10000);
                }
                else {
                    _f__SetNativityYearAndOptionsObject(value.text);
                    if (!_f__TheCurrentResearchOptionsMatchedTheStoredOnes(optionsArray)) {
                        /*  _f__infoP("Used generated theme data for "+currentRecordName); */
                        /* 	//debugP("Used generated Path!"); */
                        getThemeValues(positions[0], positions[1], positions[2], positions[3], positions[4], positions[5], positions[6], positions[7], positions[8], positions[9], positions[10], positions[11]);
                        LocalTheme = theme.slice();
                        _f__CreateOptionsObjectWithPerRecordDataIntact(optionsArray);
                    }
                    var sArrayTcopy = [
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], SubjectProcessedData
                    ];
                    var tArrayTcopy = [
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], LocalTheme
                    ];
                    if (method < 3) // 2 and 3 are external to engine and control graph display
                    {
                        var profileData = xProfile(sArrayTcopy, tArrayTcopy);
                        var barColour = profileData.pCoincidence;
                        if (!subjectFlag) {
                            inflectionPoints.push(profileData.sInflections)
                            inflectionPoints.push(profileData.sArrayTp);
                            inflectionPoints.push(profileData.sArray);
                        }
                        inflectionPoints.push(profileData.tInflections)
                        inflectionPoints.push(profileData.tArrayTp);
                        inflectionPoints.push(profileData.tArray);
                        subjectFlag = 1; // hack to ensure subject only pushed once				
                        var numberToSave = profileData.fitC + profileData.fitS; // this is the actual 'fit'
                        if ((profileData.fitC == 1) && (profileData.fitS == -1)) {
                            debugP("Error: " + currentRecordName + " no corresponding theme values outside +/- std. dev.");
                            debugP("See type 1 charts: unable to calculate fit value");
                        } else { // do not put error case names on list, pass btAccurate for graph colour
                            profileValues.push([currentRecordName, numberToSave, barColour, CurrentTimeOfBirthAccurate.toString()]);
                        }
                        subReporter = _f__pad(currentRecordName, 30, " ") + _f__pad(numberToSave.toPrecision(SAVED_NUMBER_PRECISION), 10, " ");
                        reporter[reporter.length] = subReporter;
                    }
                }
            }
        });
        if (totalFail > 0) {
            return;
        }
        profileValues.sort(function (a, b) {
            return b[1] - a[1];
        });
        for (n = 0; n < profileValues.length; n++)
            fitArray.push(profileValues[n][1]);
        if (method < 3) { // allow meaninglessly small data sets e.g. 1 item, no warning
            var fitStats = {
                mean: 0,
                sd: 0,
                sError: 0,
                skew: 0
            };
            statAnalyse(fitArray, 0, fitArray.length, fitStats, 1); // find sample sd not pop sd
            var IsInLimits = 0;
            var mean = fitStats.mean;
            var compareValueLow = fitStats.mean - (fitStats.sd * sigmaFactor);
            var compareValueHigh = fitStats.mean + (fitStats.sd * sigmaFactor);
            var TotalTest = 0;
            var OutputString = "";
            RecordToSort.header = '\n' + "Test Subject: " + $('#SubjectName').val() + ' / ' + ListToCall   +'\n'+" mean: " + pNum(fitStats.mean, precision) + " lo: " + pNum(compareValueLow, precision) + " hi: " + pNum(compareValueHigh, precision) + '\n' + '\n';
            debugP(ListToCall   +' list: '+" mean: " + pNum(fitStats.mean, precision) + " lo: " + pNum(compareValueLow, precision) + " hi: " + pNum(compareValueHigh, precision));
            var didMatchSymbol = '+';
            var didntMatchSymbol = '+'; // redefined by helper below to be actually visually useful.
            for (var i = 0; i < storeFirstData.length; i++) {
              var didMatch = false;
                 if ((storeFirstData[i][1] >= compareValueLow) && (storeFirstData[i][1] <= compareValueHigh)) {
                     didMatch = true;
                     IsInLimits++;
                     TotalIsInLimits++;
                 } else if (storeFirstData[i][1] < compareValueLow) {}
                TotalTest++;
                TotalTotalTests++;
                OutputString += "MATCH: " + didMatch + " Name: " + storeFirstData[i][0] + '&emsp;' + "mean: " + pNum(storeFirstData[i][1], precision) + '\n';
                RecordToSort.sortValuesWithString[i] = [storeFirstData[i][1], didMatch + " " + IsInLimits + "/" + TotalTest + ' ' + storeFirstData[i][0] + " fit: " + pNum(storeFirstData[i][1], precision) + '\n'];
            }
            Records[Records.length] = RecordToSort;
            TotalOutputString += OutputString;
            debugP("TotalTest " + TotalTest + " Is In Limits: " + IsInLimits);
            debugP("----------------------");
        }
       
        LoadChosenSelection(StoredSelection);
    }





function _f__ClickedxProfileLean(ListToCall) {
debugP("I CALLED _f__ClickedxProfileLean (!) method"+method);
		
		
		var reservedPrintout=[];
		
        var RecordToSort = {
            'sortValuesWithString': [],
            'header': ""
        };
        var StoredSelection = $("#selection_name").val();
        if (!($('nav select option:selected').length))
            return;
       
        var ChosenFilename = "PT_" + $('#SubjectName').val() + "_" + $("#selection_name").val() + "_" + timeStamp();
       
        if (!($('nav select option:selected').length))
            return;
        var subjectFlag = 0; // global hack to indicate 1st record in list (the subject)
        if (chartType == 1) // display purposes only, type 0 needed here
        {
            _f__warnP("select chart type 0 in Options (ALT+SHIFT+o)");
         
            return;
        }
        if (typeof (SubjectProcessedData) === "undefined" || !_f__RecordEntriesAreNoneZero(SubjectsValues)) {
            _f__warnP("Load a SUBJECT to profile");
            return;
        }
        recordsList = document.getElementById('person');
        var index = $('#person option:selected');
        if (index.length < 1) {
            _f__warnClear();
            _f__warnP("Please select TARGET record(s) from the list to cross-profile");
            return;
        }

        _f__CreateOptionsObject();
        var fitArray = [];
        var profileValues = [];
        _f__DoCompleteOptionsSaveThenLoadCycle();
        var selectedIndexLength = index.length;
        var currentSubjectName = $('#SubjectName').val().toString();
        var currentRecordName = "";
        var reprocessedRecords = [];
        var totalFail = 0;
        var br = "<br />";
        var subReporter = "Chart for " + currentSubjectName + CRLF;
        var reporter = [];
        reporter[reporter.length] = subReporter;
        subReporter = "";


        index.each(function (innerIndex, value) {
            var $this = $(this);
            currentRecordName = value.text;
            if (currentSubjectName == currentRecordName) {
                if (selectedIndexLength < 2) {
                    totalFail += 1;
                }
                return;
            }




            _f__LoadSubjectChartData(currentSubjectName);
            if ($this.length) {
                var positions = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                var LocalTheme = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                var optionsArray = new Array(OPTIONS_ARRAY_SIZE); // num options?
                var CurrentTimeOfBirthAccurate = false;
                _f__StorageRetrieveOperation(value.text, positions, LocalTheme, optionsArray);
                CurrentTimeOfBirthAccurate = _f__BirthTimeAccurate(value.text);
                if ($('#TimeImportant').prop("checked") && CurrentTimeOfBirthAccurate == false) { //debugP("Skipping record " + value.text + " for inaccurate TOB");
                    return;
                }
                if (!_f__RecordEntriesAreNoneZero(positions)) {
                    _f__infoP("I Booted the record " + value.text + " It has no meaningful POSITIONS data ", 10000);
                }
                else {
                    _f__SetNativityYearAndOptionsObject(value.text);
                    if (!_f__TheCurrentResearchOptionsMatchedTheStoredOnes(optionsArray)) {
                        /*  _f__infoP("Used generated theme data for "+currentRecordName); */
                        /* 	//debugP("Used generated Path!"); */
                        getThemeValues(positions[0], positions[1], positions[2], positions[3], positions[4], positions[5], positions[6], positions[7], positions[8], positions[9], positions[10], positions[11]);
                        LocalTheme = theme.slice();
                        _f__CreateOptionsObjectWithPerRecordDataIntact(optionsArray);
                    }
                    var sArrayTcopy = [
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], SubjectProcessedData
                    ];
                    var tArrayTcopy = [
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], LocalTheme
                    ];
                    if (method < 3) // 2 and 3 are external to engine and control graph display
                    {
                        var profileData = xProfile(sArrayTcopy, tArrayTcopy);
                        var barColour = profileData.pCoincidence;
                        if (!subjectFlag) {
                            inflectionPoints.push(profileData.sInflections)
                            inflectionPoints.push(profileData.sArrayTp);
                            inflectionPoints.push(profileData.sArray);
                        }
                        inflectionPoints.push(profileData.tInflections)
                        inflectionPoints.push(profileData.tArrayTp);
                        inflectionPoints.push(profileData.tArray);
                        subjectFlag = 1; // hack to ensure subject only pushed once				
                        var numberToSave = profileData.fitC + profileData.fitS; // this is the actual 'fit'
                        if ((profileData.fitC == 1) && (profileData.fitS == -1)) {
                            debugP("Error: " + currentRecordName + " no corresponding theme values outside +/- std. dev.");
                            debugP("See type 1 charts: unable to calculate fit value");
                        } else { // do not put error case names on list, pass btAccurate for graph colour
                            profileValues.push([currentRecordName, numberToSave, barColour, CurrentTimeOfBirthAccurate.toString()]);
                        }
                        subReporter = _f__pad(currentRecordName, 30, " ") + _f__pad(numberToSave.toPrecision(SAVED_NUMBER_PRECISION), 10, " ");
                        reporter[reporter.length] = subReporter;
                    }
					
					
                }
            }
        });
        if (totalFail > 0) {
            return;
        }
        profileValues.sort(function (a, b) {
            return b[1] - a[1];
        });
        for (n = 0; n < profileValues.length; n++)
            fitArray.push(profileValues[n][1]);
        if (method < 3) { // allow meaninglessly small data sets e.g. 1 item, no warning
            var storeFirstData = [];
            for (var i = 0; i < profileValues.length; i++) {
                storeFirstData[storeFirstData.length] = [profileValues[i][0], profileValues[i][1]];
            }
        }
        LoadChosenSelection(ListToCall);
        var subjectFlag = 0; // global hack to indicate 1st record in list (the subject)
        recordsList = document.getElementById('person');
        var index = $('#person option:selected');
        _f__CreateOptionsObject();
        var fitArray = [];
        var profileValues = [];
       
        _f__DoCompleteOptionsSaveThenLoadCycle();
        var selectedIndexLength = index.length;
        var currentSubjectName = $('#SubjectName').val().toString();
        var currentRecordName = "";
        var reprocessedRecords = [];
        var totalFail = 0;
        var br = "<br />";
        var subReporter = "Chart for " + currentSubjectName + CRLF;
        var reporter = [];
        reporter[reporter.length] = subReporter;
        subReporter = "";
        index.each(function (innerIndex, value) {
            var $this = $(this);
            currentRecordName = value.text;
            _f__LoadSubjectChartData(currentSubjectName);
            if ($this.length) {
                var positions = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                var LocalTheme = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                var optionsArray = new Array(OPTIONS_ARRAY_SIZE); // num options?
                var CurrentTimeOfBirthAccurate = false;
                _f__StorageRetrieveOperation(value.text, positions, LocalTheme, optionsArray);
                CurrentTimeOfBirthAccurate = _f__BirthTimeAccurate(value.text);
                if ($('#TimeImportant').prop("checked") && CurrentTimeOfBirthAccurate == false) {
                    return;
                }
                if (!_f__RecordEntriesAreNoneZero(positions)) {
                    _f__infoP("I Booted the record " + value.text + " It has no meaningful POSITIONS data ", 10000);
                }
                else {
                    _f__SetNativityYearAndOptionsObject(value.text);
                    if (!_f__TheCurrentResearchOptionsMatchedTheStoredOnes(optionsArray)) {
                        /*  _f__infoP("Used generated theme data for "+currentRecordName); */
                        /* 	//debugP("Used generated Path!"); */
                        getThemeValues(positions[0], positions[1], positions[2], positions[3], positions[4], positions[5], positions[6], positions[7], positions[8], positions[9], positions[10], positions[11]);
                        LocalTheme = theme.slice();
                        _f__CreateOptionsObjectWithPerRecordDataIntact(optionsArray);
                    }
                    var sArrayTcopy = [
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], SubjectProcessedData
                    ];
                    var tArrayTcopy = [
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], LocalTheme
                    ];
                    if (method < 3) // 2 and 3 are external to engine and control graph display
                    {
                        var profileData = xProfile(sArrayTcopy, tArrayTcopy);
                        var barColour = profileData.pCoincidence;
                        if (!subjectFlag) {
                            inflectionPoints.push(profileData.sInflections)
                            inflectionPoints.push(profileData.sArrayTp);
                            inflectionPoints.push(profileData.sArray);
                        }
                        inflectionPoints.push(profileData.tInflections)
                        inflectionPoints.push(profileData.tArrayTp);
                        inflectionPoints.push(profileData.tArray);
                        subjectFlag = 1; // hack to ensure subject only pushed once				
                        var numberToSave = profileData.fitC + profileData.fitS; // this is the actual 'fit'
                        if ((profileData.fitC == 1) && (profileData.fitS == -1)) {
                            debugP("Error: " + currentRecordName + " no corresponding theme values outside +/- std. dev.");
                            debugP("See type 1 charts: unable to calculate fit value");
                        } else { // do not put error case names on list, pass btAccurate for graph colour
                            profileValues.push([currentRecordName, numberToSave, barColour, CurrentTimeOfBirthAccurate.toString()]);
                        }
                        subReporter = _f__pad(currentRecordName, 30, " ") + _f__pad(numberToSave.toPrecision(SAVED_NUMBER_PRECISION), 10, " ");
                        reporter[reporter.length] = subReporter;
                    }
                }
            }
        });


        if (totalFail > 0) {
            return;
        }
        profileValues.sort(function (a, b) {
            return b[1] - a[1];
        });
        for (n = 0; n < profileValues.length; n++)
            fitArray.push(profileValues[n][1]);
        if (method < 3) { // allow meaninglessly small data sets e.g. 1 item, no warning
            var fitStats = {
                mean: 0,
                sd: 0,
                sError: 0,
                skew: 0
            };
            statAnalyse(fitArray, 0, fitArray.length, fitStats, 1); // find sample sd not pop sd
            var IsInLimits = 0;
            var mean = fitStats.mean;
var compareValueLow = fitStats.mean - fitStats.sd;
var compareValueHigh = fitStats.mean + fitStats.sd;
            var TotalTest = 0;
            var OutputString = "";
            RecordToSort.header = '\n' + "Test Subject: " + $('#SubjectName').val() + ' / ' + ListToCall   +'\n'+" mean: " + pNum(fitStats.mean, precision) + " lo: " + pNum(compareValueLow, precision) + " hi: " + pNum(compareValueHigh, precision) + '\n' + '\n';
            debugP(ListToCall   +' list: '+" mean: " + pNum(fitStats.mean, precision) + " lo: " + pNum(compareValueLow, precision) + " hi: " + pNum(compareValueHigh, precision));
            var didMatchSymbol = '+';
            var didntMatchSymbol = '+'; // redefined by helper below to be actually visually useful.
            for (var i = 0; i < storeFirstData.length; i++) {
              var didMatch = false;
                 if ((storeFirstData[i][1] >= compareValueLow) && (storeFirstData[i][1] <= compareValueHigh)) {
                     didMatch = true;
                     IsInLimits++;
                     TotalIsInLimits++;
                 } else if (storeFirstData[i][1] < compareValueLow) {}
                TotalTest++;
                TotalTotalTests++;
                debugP(didMatch + " " + IsInLimits + "/" + TotalTest + '&emsp;' + "fit: " + pNum(storeFirstData[i][1], precision) + "&emsp;" + storeFirstData[i][0]); // + " lo: "+pNum(compareValueLow, precision) + " hi: "+pNum(compareValueHigh,precision)    );
                OutputString += "MATCH: " + didMatch + " Name: " + storeFirstData[i][0] + '&emsp;' + "mean: " + pNum(storeFirstData[i][1], precision) + '\n';
                RecordToSort.sortValuesWithString[i] = [storeFirstData[i][1], didMatch + " " + IsInLimits + "/" + TotalTest + ' ' + storeFirstData[i][0] + " fit: " + pNum(storeFirstData[i][1], precision) + '\n'];
            }
            Records[Records.length] = RecordToSort;
            TotalOutputString += OutputString;
            debugP("TotalTest " + TotalTest + " Is In Limits: " + IsInLimits);
            debugP("----------------------");
        }
       
        LoadChosenSelection(StoredSelection);
    }






function _f__ResetTheEngine() {
    Records = [];
    fitArray = [];
    SubjectValues = [];
    TotalIsInLimits = 0;
    TotalTotalTests = 0;
    TotalOutputString = "";
}

function CreateArrayCopyToPrecision(arrayToCopy)
{
	var newArray = [];
	
	for(var i =0; i < arrayToCopy.length;     i++) {
		
	factor = Math.pow ( 10, precision );
		var num = arrayToCopy[i] * factor;
		var resultant = Math.ceil ( num / (Math.pow ( 10, precision )));

	
           newArray.push( resultant    );
    
	}
	
	return newArray;
	
}


function arraysEqual(arr1, arr2) {
    if(arr1.length !== arr2.length)
        return false;
    for(var i = arr1.length; i--;) {
        if(arr1[i] !== arr2[i])
            return false;
    }

    return true;
}

function _f__RetrieveThemeValueForNameAndCompare(name, compareValues)
{
	var storedValues = _f__decode(localStorage.getItem(name)).split(delimeter);
	 var  retrievedThemes =[];
	var simpleIndex = 0;
        for (i = 0; i < 24; i += 2) {
          
            retrievedThemes[simpleIndex] = parseFloat(storedValues[i + 1]);
            simpleIndex++;
        }
	
	
	var test1 = CreateArrayCopyToPrecision(retrievedThemes);
	var test2 = CreateArrayCopyToPrecision(compareValues);
	
	return arraysEqual(test1,  test2);
	
	
}




var TotalIsInLimits = 0;
var TotalTotalTests = 0;
var TotalOutputString = "";
function _f__ClickedxProfile() {
debugP("I CALLED _f__ClickedxProfile() method"+method);
    if (method == 3) {
        method = 1; // set default method since this option does not call profiling functions but
        if (!(inflectionPoints.length))
            return;
        if (!($('nav select option:selected').length))
            return;
        var ChosenFilename = "";
        if (!($("#selections li").hasClass("selected"))) {
            var chosenInput = prompt("Enter filename");
            ChosenFilename = chosenInput;
            SaveCurrentSelectionWithName(ChosenFilename);
        }
        if (ChosenFilename == "")
            ChosenFilename = $("#selection_name").val();
        if (!($('nav select option:selected').length))
            return;
        curveDataSet = fnCoeffs(inflectionPoints);
        var rr = [];
        rr[0] = $('#SubjectName').val();
        $('nav select option:selected').each(function (i, selected) {
            rr[i + 1] = $(selected).text();
        });
        _f__createDownloadableFile(ChosenFilename + "_names.txt", JSON.stringify(rr.join('\n')));
        _f__createDownloadableFile(ChosenFilename + "_data.txt", JSON.stringify(curveDataSet));
        return;
    }
    var subjectFlag = 0; // global hack to indicate 1st record in list (the subject)
    if (chartType == 1) // display purposes only, type 0 needed here
    {
        _f__warnP("select chart type 0 in Options (ALT+SHIFT+o)");
        return;
    }
    if (typeof (SubjectProcessedData) === "undefined" || !_f__RecordEntriesAreNoneZero(SubjectsValues)) {
        _f__warnP("Load a SUBJECT to profile");
        return;
    }
    recordsList = document.getElementById('person');
    var index = $('#person option:selected');
    if (index.length < 1) {
        _f__warnClear();
        _f__warnP("Please select TARGET record(s) from the list to cross-profile");
        return;
    }
    _f__CreateOptionsObject();
    var fitArray = [];
    var profileValues = [];
    _f__clearCharts();
    _f__debugClear(); // Do NOT clear debug proactively. WJ18 request 31/01/06 (reinstted 22/03/2018 - Will)
    _f__warnClear();
    _f__DoCompleteOptionsSaveThenLoadCycle();
    var selectedIndexLength = index.length;
    var currentSubjectName = $('#SubjectName').val().toString();
    var currentRecordName = "";
    var reprocessedRecords = [];
    var totalFail = 0;
    var br = "<br />";
    var subReporter = "Chart for " + currentSubjectName + CRLF;
    var reporter = [];
    reporter[reporter.length] = subReporter;
    subReporter = "";
    index.each(function (innerIndex, value) {
        var $this = $(this);
        currentRecordName = value.text;
        if (currentSubjectName == currentRecordName) {
            if (selectedIndexLength < 2) {
                totalFail += 1;
                _f__warnP("DONT select the SUBJECT as a TARGET from the list to cross-profile");
            }
            return;
        }
        _f__LoadSubjectChartData(currentSubjectName);
        if ($this.length) {
            var positions = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            var LocalTheme = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            var optionsArray = new Array(OPTIONS_ARRAY_SIZE); // num options?
            var CurrentTimeOfBirthAccurate = false;
            _f__StorageRetrieveOperation(value.text, positions, LocalTheme, optionsArray);
            CurrentTimeOfBirthAccurate = _f__BirthTimeAccurate(value.text);
            if ($('#TimeImportant').prop("checked") && CurrentTimeOfBirthAccurate == false) { //debugP("Skipping record " + value.text + " for inaccurate TOB");
                return;
            }
            if (!_f__RecordEntriesAreNoneZero(positions)) {
                _f__infoP("I Booted the record " + value.text + " It has no meaningful POSITIONS data ", 10000);
            }
            else {
                _f__SetNativityYearAndOptionsObject(value.text);
                if (!_f__TheCurrentResearchOptionsMatchedTheStoredOnes(optionsArray)) {
                    /*  _f__infoP("Used generated theme data for "+currentRecordName); */
                    /* 	//debugP("Used generated Path!"); */
              

			


			
					

					getThemeValues(positions[0], positions[1], positions[2], positions[3], positions[4], positions[5], positions[6], positions[7], positions[8], positions[9], positions[10], positions[11]);
                    LocalTheme = theme.slice();
					
					
					
					
                    /*
                    debugP(LocalTheme);
                    debugP("record.js: LocalTheme");
                    for (n=0; n<12; n++)
                    	debugP(pNum(LocalTheme[n], precision));
                    */
                    _f__CreateOptionsObjectWithPerRecordDataIntact(optionsArray);
                    /*
                                    if (!localStorage.getItem(currentRecordName)) {
                                        _f__StorageSaveOperation(currentRecordName, positions, LocalTheme, GlobalCurrentSaveOptions);
                                    } else {
                                        _f__deleteRecordFromRecordsList(currentRecordName);
                                        _f__StorageSaveOperation(currentRecordName, positions, LocalTheme, GlobalCurrentSaveOptions);
                                    }
                    */
                }
                var sArrayTcopy = [
                    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], SubjectProcessedData
                ];
                var tArrayTcopy = [
                    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], LocalTheme
                ];

				
				
				
				
				
				
				/*
                debugP("record.js: sArrayCopy");
                for (n=0; n<12; n++)
                	debugP(pNum(sArrayTcopy[1][n], precision));
                debugP("record.js: tArrayCopy");
                for (n=0; n<12; n++)
                	debugP(pNum(tArrayTcopy[1][n], precision));
*/              
                if (method < 3) // 2 and 3 are external to engine and control graph display
                {
                    var profileData = xProfile(sArrayTcopy, tArrayTcopy);
                    var barColour = profileData.pCoincidence;
                    if (!subjectFlag) {
                        inflectionPoints.push(profileData.sInflections)
                        inflectionPoints.push(profileData.sArrayTp);
                        inflectionPoints.push(profileData.sArray);
                    }
                    inflectionPoints.push(profileData.tInflections)
                    inflectionPoints.push(profileData.tArrayTp);
                    inflectionPoints.push(profileData.tArray);
                    subjectFlag = 1; // hack to ensure subject only pushed once				
                    var numberToSave = profileData.fitC + profileData.fitS; // this is the actual 'fit'
                    if ((profileData.fitC == 1) && (profileData.fitS == -1)) {
                        debugP("Error: " + currentRecordName + " no corresponding theme values outside +/- std. dev.");
                        debugP("See type 1 charts: unable to calculate fit value");
                    } else { // do not put error case names on list, pass btAccurate for graph colour
                        profileValues.push([currentRecordName, numberToSave, barColour, CurrentTimeOfBirthAccurate.toString()]);
                    }
                    subReporter = _f__pad(currentRecordName, 30, " ") + _f__pad(numberToSave.toPrecision(SAVED_NUMBER_PRECISION), 10, " ");
                    reporter[reporter.length] = subReporter;
                }
				
				
					
				
            }
			
			
				
					
			
        }
    });
    if (totalFail > 0) {
        return;
    }
    profileValues.sort(function (a, b) {
        return b[1] - a[1];
    });
    for (n = 0; n < profileValues.length; n++)
        fitArray.push(profileValues[n][1]);
    if (method < 3) { // allow meaninglessly small data sets e.g. 1 item, no warning
        var fitStats = {
            mean: 0,
            sd: 0,
            sError: 0,
            skew: 0
        };
        statAnalyse(fitArray, 0, fitArray.length, fitStats, 1); // find sample sd not pop sd

var count;
for ( count = 0; count < fitArray.length; count++ )

if  (count > 1){
	profileValues[count][1] -= fitStats.mean;	// force zero offset
}

        _f__drawxProfileChart("xProfile chart", profileValues, 'chart_div_4', fitStats);
    }
    /* var br = "<br />"; */
    /* reporter += profileValues; */
    lastFullChartReport = reporter;
    /* 
	var realsdFactor =$(window).sdFactor();
	var modulus = 1;
	
	if( realsdFactor < 1080)
       modulus = 1;

    if (realsdFactor < 1500)
        modulus = 2;

	if( realsdFactor < 1890)
        modulus = 3;
	if ( realsdFactor > 1890)
        modulus = 4;
	
    for (l = 0; l < profileValues.length; l++) {

		var temp = profileValues[l][0] + "&nbsp;" + profileValues[l][1];

	reporter +=	_f__pad( temp, 57, "&nbsp;");	
		
		if (l % modulus ==0 )
		{  
	reporter += br;
	}
	else
	{
		reporter+= ", ";
	}
    } */
    /* //debugP(reporter); */
    /* 	lastFullChartReport = profileValues; */
	
	_f__ResetTheEngine();
	
}

function _f__CreateOptionsObjectWithPerRecordDataIntact(optionsArray) {
        var passeddatetime = optionsArray[2];
        var passedtimeAccurate = optionsArray[3];
        var passedlatitude = optionsArray[1];
        var passedlongitude = optionsArray[0];
        var passedUt = optionsArray[4];
        var passedDst = optionsArray[5];
        var passedaoIndex = optionsArray[6];
        var passedtfIndex = optionsArray[7];
        var passedprecessionFlag = optionsArray[8];
        var passedorbType = optionsArray[9];
        var passedpoindex = optionsArray[10];
        var passedcharttype = optionsArray[11];
        _f__CreateOptionsObjectButPreservePerUserRecordVariablesToCurrentState();
        GlobalCurrentSaveOptions[0] = passedlongitude;
        GlobalCurrentSaveOptions[1] = passedlatitude;
        GlobalCurrentSaveOptions[2] = passeddatetime; // DATETIME of BIRTH
        GlobalCurrentSaveOptions[3] = passedtimeAccurate
        GlobalCurrentSaveOptions[4] = passedUt; // seconds
        GlobalCurrentSaveOptions[5] = passedDst; // seconds
    }
function _f__populateHMSArraysFromValueArraysValues() {
    $('.numVals').each(function () {
        _f__ConvertDecimalToHMSSingle($(this));
    })
}

function _f__pad(n, sdFactor, z) {
    z = z || '0';
    n = n + '';
    return n.length >= sdFactor ? n : new Array(sdFactor - n.length + 1).join(z) + n;
}

function _f__padLeft(n, sdFactor, z) {
    z = z || '0';
    n = n + '';
    return n.length >= sdFactor ? n : n + new Array(sdFactor - n.length + 1).join(z);
}

function _f__SetNativityYearAndOptionsObject(NameOfRecord) {
        if (0 | OptionsObject.precessionFlag) {
            if (_f__GetNativityYear(NameOfRecord, nativityYear)) {
                precessionFlag = 1;
                /*  //debugP(NameOfRecord +" "+ nativityYear); */
            }
            else {
                /*	// doesn't exist in options
                	OptionsObject.precessionFlag =0;	
                	precessionFlag = 0;
                */
                /* //debugP(NameOfRecord +" "+ nativityYear); */
            }
        }
    }
function _f__GetNativityYear(NameOfRecord, _s__OUT_forDate) {
        var storedValues = [];
        var dateIndexInStoredRecords = 0;
        if (localStorage.getItem(NameOfRecord) === null) {
            _s__OUT_forDate = nativityYear = parseInt('0000');
            return false;
        }
        else {
            var storedValues = _f__decode(localStorage.getItem(NameOfRecord)).split(delimeter);
            var dateIndexInStoredRecords = 26;
            var splitDate = storedValues[dateIndexInStoredRecords].trim().replace(/\//g, ' ').replace(/:/g, ' ').split(" ");
            _s__OUT_forDate = nativityYear = parseInt(splitDate[2]);
        }
        var fullDOBcheck = splitDate[0] + ' ' + splitDate[1] + ' ' + splitDate[2];
        if (_s__OUT_forDate == '0000' && fullDOBcheck == __DOB_DEFAULTS__) {
            return false;
        }
        return true;
    }
function _f__SetNativityYear(NameOfRecord) {
        var storedValues = [];
        var dateIndexInStoredRecords = 0;
        if (localStorage.getItem(NameOfRecord) === null) {
            var splitDate = OptionsObjectDateTime.trim().replace(/\//g, ' ').replace(/:/g, ' ').split(" ");
            nativityYear = parseInt(splitDate[0]);
        }
        else {
            var storedValues = _f__decode(localStorage.getItem(NameOfRecord)).split(delimeter);
            var dateIndexInStoredRecords = 26;
            var splitDate = storedValues[dateIndexInStoredRecords].trim().replace(/\//g, ' ').replace(/:/g, ' ').split(" ");
            nativityYear = parseInt(splitDate[2]);
        }
    }
function _f__populateValueArraysFromStoredValues(NameOfRecord, NameOfColumn) {
        var IsSubjectSave = false;
        if (NameOfColumn == "subject") {
            IsSubjectSave = true;
        }
        toastr.clear();
        var MangledDataDetectedTidyUpRequired = false;
        var storedValues = _f__decode(localStorage.getItem(NameOfRecord)).split(delimeter);
        _f__InitProcessedDataArrayToBlanksIfUndefined();
        var simpleIndex = 0;
        var i = 24;
        OptionsObjectLon = parseFloat(storedValues[i]);
        if (parseFloat(OptionsObjectLon).IsDefault() ||
            OptionsObjectLon.toString().IsDefault()) {
            OptionsObjectLon = __LON_DEFAULTS__.toString();
            MangledDataDetectedTidyUpRequired = true;
            /*   //debugP("Fixed Mangled Lon");
	alert(OptionsObjectLon); */
        }
        storedValues[i] = OptionsObjectLon;
        lat = storedValues[i];
        /* alert(lat); */
        if (IsSubjectSave) {
            $('#uk-lon').val(OptionsObjectLon);
        } else {
            $('#uk-lon-tar').val(OptionsObjectLon);
        }
        i++;
        OptionsObjectLat = parseFloat(storedValues[i]);
        if (parseFloat(OptionsObjectLat).IsDefault() ||
            OptionsObjectLat.toString().IsDefault()) {
            OptionsObjectLat = __LAT_DEFAULTS__.toString();
            /*   //debugP("Fixed Mangled Lat"); */
            MangledDataDetectedTidyUpRequired = true;
        }
        storedValues[i] = OptionsObjectLat;
        lon = storedValues[i];
        /* 	alert(lon);	 */
        if (IsSubjectSave) {
            $('#uk-lat').val(OptionsObjectLat);
        } else {
            $('#uk-lat-tar').val(OptionsObjectLat);
        }
        i++;
        var splitDate = storedValues[i].trim().replace(/\//g, ' ').replace(/:/g, ' ').split(" ");
        var activeDateNode = "";
        if (IsSubjectSave) {
            activeDateNode = $('#uk-dob');
        } else {
            activeDateNode = $('#uk-dob-tar');
        }
        if (typeof splitDate[1] === "undefined" || typeof splitDate[2] === "undefined") {
            /* 	//debugP("set dob defs");  */
            MangledDataDetectedTidyUpRequired = true; // see, told ya.  Set the flag to save back the record.
            splitDate[1] = "00"; // NB. preceding a number with a zero would make it octal.. dont do that.
            splitDate[2] = "0000";
            activeDateNode.val(__DOB_DEFAULTS__);
            dob = __DOB_DEFAULTS__;
        } else {
            var tidiedVal = splitDate[0] + ' ' + splitDate[1] + ' ' + splitDate[2];
            activeDateNode.val(tidiedVal);
            dob = tidiedVal;
        }
        var activeTimeNode = "";
        if (IsSubjectSave) {
            activeTimeNode = $('#uk-tob');
        } else {
            activeTimeNode = $('#uk-tob-tar');
        }
        if (typeof splitDate[3] === "undefined" || typeof splitDate[4] === "undefined") {
            /* //debugP("set tob defs"); */
            MangledDataDetectedTidyUpRequired = true; // see, told ya.  Set the flag to save back the record.
            splitDate[3] = "00";
            splitDate[4] = "00";
            activeTimeNode.val(__TOB_DEFAULTS__);
            tob = __TOB_DEFAULTS__;
        } else {
            var tidiedVal = splitDate[3] + ' ' + splitDate[4];
            activeTimeNode.val(tidiedVal);
            tob = tidiedVal;
        }
        OptionsObjectDateTime = new Date(splitDate[2], splitDate[1], splitDate[0], splitDate[3], splitDate[4]);
        var datetime = OptionsObjectDateTime;
        storedValues[i] =
            splitDate[0] + " " +
            splitDate[1] + " " +
            splitDate[2] + " " +
            splitDate[3] + " " +
            splitDate[4];
        i++;
        var activeAccNode = "";
        if (IsSubjectSave) {
            activeAccNode = $('#uk-acc')
        } else {
            activeAccNode = $('#uk-acc-tar');
        }
        OptionsObjectDateTimeAccurate = parseInt(storedValues[i]);
        var isAccurate = OptionsObjectDateTimeAccurate;
        activeAccNode.prop("checked", isAccurate);
        storedValues[i] = isAccurate;
        acc = isAccurate;
        i++;
        OptionsObjectUtOffset = parseInt(storedValues[i]);
        i++;
        OptionsObjectDstOffset = parseInt(storedValues[i]);
        i++;
        $.each(OptionsObject, function (index, value) {
            OptionsObject[index] = storedValues[i];
            i++;
            /*   	$('#$(index)').val(storedValues[i]); */
        });
        _f__AppendPreparedStateNotSetFromWidgetsThenCreateOptions(NameOfColumn);
        var temporaryCopyOfGlobalCurrentSaveOptions = GlobalCurrentSaveOptions.slice();
        var itDoesntMatch = !_f__TheCurrentResearchOptionsMatchedTheStoredOnes(temporaryCopyOfGlobalCurrentSaveOptions);
        var itDoesMatch = !itDoesntMatch;
        _f__SetNativityYearAndOptionsObject(NameOfRecord);
        /* 	alert(storedValues); */
        if (itDoesntMatch || MangledDataDetectedTidyUpRequired) {
            /* //debugP("no match"); */
            if (NameOfColumn == "subject") {
                _f__LoadSubjectChartData(NameOfRecord);
            } else {
                _f__LoadTargetChartData(NameOfRecord);
            }
            for (internal = 0; internal < 24; internal += 2) {
                if (NameOfColumn == "subject") {
                    SubjectsValues[simpleIndex] = document.f1.numvalue1[simpleIndex].value = parseFloat(storedValues[internal]).toFixed(5);
                    $("#SubjectName").val(NameOfRecord);
                } else {
                    TargetsValues[simpleIndex] = document.f1.numvalue2[simpleIndex].value = parseFloat(storedValues[internal]).toFixed(5);
                    $("#TargetName").val(NameOfRecord);
                }
                simpleIndex++;
            }
            _f__SetOptionsObjectFromWidgetState(NameOfColumn);
            var fixedAmount = 6;
            /* alert(temporaryCopyOfGlobalCurrentSaveOptions); */
            GlobalCurrentSaveOptions = temporaryCopyOfGlobalCurrentSaveOptions.slice(0, fixedAmount);
            /* 		//debugP("PreGLobal: "+GlobalCurrentSaveOptions); */
            GlobalCurrentSaveOptions[0] = OptionsObjectLon = storedValues[24];
            GlobalCurrentSaveOptions[1] = OptionsObjectLat = storedValues[25];
            /* 		//debugP("PostGLobal: "+GlobalCurrentSaveOptions); */
            _f__AppendResearchFlagsToGlobalOptionsObject();
            /* 		//debugP("PostGLobalAppend: "+GlobalCurrentSaveOptions); */
            /* 	_f__infoP('Updated record for '+NameOfRecord+' replacing record:\n'+temporaryCopyOfGlobalCurrentSaveOptions+'<br/> With '+GlobalCurrentSaveOptions+' for column '+NameOfColumn);  */
            _f__doConditionalSave(NameOfColumn, NameOfRecord, GlobalCurrentSaveOptions);
        } else {
            for (internal = 0; internal < 24; internal += 2) {
                if (NameOfColumn == "subject") {
                    SubjectsValues[simpleIndex] = document.f1.numvalue1[simpleIndex].value = parseFloat(storedValues[internal]).toFixed(5);
                    SubjectProcessedData[simpleIndex] = parseFloat(storedValues[internal + 1]).toFixed(4);
                    $("#SubjectName").val(NameOfRecord);
                } else {
                    TargetsValues[simpleIndex] = document.f1.numvalue2[simpleIndex].value = parseFloat(storedValues[internal]).toFixed(STORED_VALUES_INTERNAL_PRECISION);
                    TargetsProcessedData[simpleIndex] = parseFloat(storedValues[internal + 1]).toFixed(TARGET_PROCESSED_DATA_PRECISION);
                    $("#TargetName").val(NameOfRecord);
                }
                simpleIndex++;
            }
        }
        /* 	//debugP("By now we should have a stored record that reflects the current settings."); */
    }
function _f__deleteRecordFromSelectionList(name) {
        if (localStorage.getItem(name)) {
            localStorage.removeItem(name)
        };
        _f__populateSelectionsList();
    }
function _f__deleteRecordFromRecordsList(name) {
    if (localStorage.getItem(name)) {
        localStorage.removeItem(name)
    };
    _f__populateRecordsList();
}

function _f__sortLocalStorage() {
        if (localStorage.length > 0) {
            var localStorageArray = new Array();
            for (i = 0; i < localStorage.length; i++) {
                localStorageArray[i] = localStorage.key(i) + localStorage.getItem(localStorage.key(i));
            }
        }
        var sortedArray = localStorageArray.sort();
        return sortedArray;
    }
  



function _f__populateRecordsList() {
	
	
	
        recordsList = document.getElementById('person');
        var OptionsObjectPerRecordStructure = ['Long', 'Lat', 'Date/Time', 'TimeAccurate', 'Offset', 'DstOffset'];
        var OptionsObjectCompleteStructure = [];
        jQuery.each(OptionsObjectPerRecordStructure, function (key, value) {
            OptionsObjectCompleteStructure[OptionsObjectCompleteStructure.length] = value;
        });
        jQuery.each(OptionsObject, function (key, value) {
            OptionsObjectCompleteStructure[OptionsObjectCompleteStructure.length] = key.toString();
        });
        var total = 0;
        if (typeof SubjectsValues !== "undefined") {
            for (var i = 0; i < 12; i++) {
                total += SubjectsValues[i]
            }
        }
        if (total <= 0) {
            $('#SubjectName').val("Subject");
        }
		
		$(recordsList).empty();
		
        var name = new Array();
        if (_f__storageAvailable('localStorage') == false) {
            _f__warnP('Too bad, no localStorage for us. ');
            return;
        }
        if (localStorage.length == 0) {
            _f__InitProcessedDataArrayToBlanksIfUndefined();
            /*  #### ADD SOME DEFAULT STUFF FOR TESTING
                 if (!localStorage.getItem('ExampleRecord' + localStorage.length)) {
                  StorageSaveOperation('ExampleRecord' + localStorage.length, SubjectsValues, SubjectProcessedData);
                   }
               if (!localStorage.getItem('ExampleOtherRecord' + localStorage.length)) {
                    StorageSaveOperation('ExampleOtherRecord' + localStorage.length, TargetsValues, TargetsProcessedData);
              } */
        }
		
var OurCompleteStoredRecordcount =0; 
		
	
for (var i = 0; i <= localStorage.length - 1; i++) {
  	
    var key = localStorage.key(i);
            var val = _f__decode(localStorage.getItem(key));
            if (key.includes(__PLAYLIST__PREFIX__)) {
                continue;
            }
            if (key.startsWith(__PREFIX__)) {
                continue;
            }
OurCompleteStoredRecordcount++;
}
		
        for (var i = 0; i <= localStorage.length - 1; i++) {
            var key = localStorage.key(i);
            var val = _f__decode(localStorage.getItem(key));
            if (key.includes(__PLAYLIST__PREFIX__)) {
                continue;
            }
            if (key.startsWith(__PREFIX__)) {
                continue;
            }
			
		

		
			
            var value = val.split(delimeter); //splitting string inside array to get name
            var opt = document.createElement('option');
            opt.value = value;
            var recordsCommaDelimited = _f__decode(localStorage.getItem(key)).split(delimeter);
            var recordsTidyForTooltips = "";
            var ind = 0;
            $.each(recordsCommaDelimited, function (key, value) {
                ++ind;
                if (key >= 24) {
                    var offsetter = (ind - 25);
                    recordsTidyForTooltips += OptionsObjectCompleteStructure[offsetter] + ": ";
                }
                var FixedLengthString = value.toString() + "                ";
                FixedLengthString = FixedLengthString.slice(0, 15);
                recordsTidyForTooltips += FixedLengthString;
                if (ind % 2 == 0 && key < 24) {
                    recordsTidyForTooltips += "\n";
                } else if (key >= 24) {
                    recordsTidyForTooltips += "\n";
                }
            })
            opt.title = recordsTidyForTooltips;
            opt.innerHTML = key;
            recordsList.appendChild(opt);
        }
        $('select.personRecords').sortSelect();
		
		
		 if (  (OurCompleteStoredRecordcount + 0 ) !=   (  $(".personRecords option").length 	 + 0)	)
		 {
			
			 _f__warnP( OurCompleteStoredRecordcount +" DOES NOT EQUAL "+   $(".personRecords option").length     );
			  debugP('test FAILED');
			 
		 }
		
    }
    /**
     * Sort values alphabetically in select
     * source: http://stackoverflow.com/questions/12073270/sorting-options-elements-alphabetically-using-jquery
     */
$.fn.extend({
    sortSelect() {
        let options = this.find("option"),
            arr = options.map(function (_, o) {
                return {
                    t: $(o).text(),
                    v: o.value,
                    x: o.title
                };
            }).get();
        arr.sort((o1, o2) => { // sort select
            let t1 = o1.t.toLowerCase(),
                t2 = o2.t.toLowerCase();
            return t1 > t2 ? 1 : t1 < t2 ? -1 : 0;
        });
        options.each((i, o) => {
            o.value = arr[i].v;
            $(o).text(arr[i].t);
            o.title = arr[i].x;
        });
    }
});
function _f__StorageSaveOperation(nameOfThing, positions, themes, options) {
    var storedValues = "";
    for (i = 0; i < 12; i++) {
        storedValues += parseFloat(positions[i]).toFixed(SAVED_NUMBER_PRECISION)
        storedValues += delimeter;
        storedValues += parseFloat(themes[i]).toFixed(SAVED_NUMBER_PRECISION);
        storedValues += delimeter;
    }
    $.each(options, function (index, value) {
        storedValues += value;
        storedValues += delimeter;
    });
    storedValues = storedValues.slice(0, -1);
    storedValues = _f__encode(storedValues);
    localStorage.setItem(nameOfThing, storedValues);
    _f__populateRecordsList();
}

function _f__BirthTimeAccurate(nameOfThing) {
        function parseBool(str) {
                if (str.length == null) {
                    return str == 1 ? true : false;
                } else {
                    return str == "1" ? true : false;
                }
            }
        var storedValues = _f__decode(localStorage.getItem(nameOfThing)).split(delimeter);

  return storedValues[27];

    }
function _f__StorageRetrieveOperation(nameOfThing, positions, themes, options) {
        TimeAccurateBool = false;
        var storedValues = _f__decode(localStorage.getItem(nameOfThing)).split(delimeter);
        var simpleIndex = 0;
        for (i = 0; i < 24; i += 2) {
          



		positions[simpleIndex] = pNum(parseFloat(storedValues[i]), precision);

     	themes[simpleIndex] =   pNum(parseFloat(storedValues[i + 1]), precision);
            simpleIndex++;
        }
        var simpleIndexForJustHere = 0;
        for (i = 24; i < storedValues.length; i++) {
            options[simpleIndexForJustHere] = storedValues[i];
            simpleIndexForJustHere++;
            simpleIndex++;
        }
    }
function _f__fnOpenDeleteDialog() {
        var index = $("#person").prop('selectedIndex');
        var chosen = $("#person :selected");
        var namer = "";
        if (index < 0) {
            _f__warnClear();
            _f__warnP("Please select a single record from the list to delete");
            return;
        }
        var texter = "Really Delete ";
        $.each(chosen, function (key, value) {
            texter += "<br />" + value.text + " ";
            namer += value.text + " ";
        })
        texter += "?";
        $("#dialog-confirm").html(texter);
        $("#dialog-confirm").show();
        $("#dialog-confirm").dialog({
            resizable: false,
            modal: true,
            title: "DELETE " + namer,
            height: 250,
            sdFactor: 400,
            buttons: {
                "No": function () {
                    $(this).dialog('close');
                    _f__warnP("record NOT deleted");
                },
                "Yes": function () {
                    $(this).dialog('close');
                    _f__ClickedRemove();
                }
            }
        }).dialog('widget').position({
            my: 'top',
            at: 'top+100',
            of: window
        });
    }
function _f__ClickedRemove() {
        var chosen = $("#person :selected");
        _f__warnClear();
        $.each(chosen, function (key, value) {
            _f__deleteRecordFromRecordsList(value.text);
        })
        _f__infoP("Deleted " + chosen.length + " records", 3000);
        var TotalContentCount = $('#person option').length;
        $('#RecordCount').text(TotalContentCount + " records")
    }
function _f__fnOpenDeleteSelectionDialog() {
    _f__warnClear();
    if (!$('#selections li').hasClass("selected")) {
        _f__warnP('Select a record to delete ');
        return;
    }
    var thingToDel = $('#selections li.selected').text();
    var namer = '';
    var texter = "Really Delete ";
    texter += "<br />" + thingToDel + " ";
    namer += thingToDel + " ";
    texter += "?";
    $("#dialog-confirm").html(texter);
    $("#dialog-confirm").show();
    $("#dialog-confirm").dialog({
        resizable: false,
        modal: true,
        title: "DELETE " + namer,
        height: 250,
        sdFactor: 400,
        buttons: {
            "No": function () {
                $(this).dialog('close');
                _f__warnP("Selection Group NOT deleted");
            },
            "Yes": function () {
                $(this).dialog('close');
                _f__ClickedRemoveSelection();
            }
        }
    }).dialog('widget').position({
        my: 'top',
        at: 'top+100',
        of: window
    });
}

function DeleteCurrentSelection() {
        var thingToDel = $('#selections li.selected').text();
        if (localStorage.getItem(__PLAYLIST__PREFIX__ + thingToDel)) {
            localStorage.removeItem(__PLAYLIST__PREFIX__ + thingToDel)
        };
        _f__populateSelectionsList();
    }
function _f__ClickedRemoveSelection() {
        DeleteCurrentSelection();
    }
function _f__RecordEntriesAreNoneZero(arrayToCheck) {
        var ret = false
        var total = 0;
        for (i = 0; i < 12; i++) {
            total += parseFloat(arrayToCheck[i])
        }
        if (total >= 1)
            ret = true;
        return ret;
    }
function _f__ClickedSave(event, NameOfColumn) {
        _f__CreateOptionsObject();
        var IsSubjectSave = false;
        if (NameOfColumn == "subject") {
            IsSubjectSave = true;
        }
        var activeLonNode = "";
        var activeLatNode = "";
        var activeDobNode = "";
        var activeTobNode = "";
        var activeAccNode = "";
        var activeNameNode = "";
        if (IsSubjectSave) {
            activeLonNode = $('#uk-lon');
            activeLatNode = $('#uk-lat');
            activeDobNode = $('#uk-dob');
            activeTobNode = $('#uk-tob');
            activeAccNode = $('#uk-acc');
            activeNameNode = $('#SubjectName');
        } else {
            activeLonNode = $('#uk-lon-tar');
            activeLatNode = $('#uk-lat-tar');
            activeDobNode = $('#uk-dob-tar');
            activeTobNode = $('#uk-tob-tar');
            activeAccNode = $('#uk-acc-tar');
            activeNameNode = $('#TargetName');
        }
        GlobalCurrentSaveOptions[0] = activeLonNode.val();
        GlobalCurrentSaveOptions[1] = activeLatNode.val();
        GlobalCurrentSaveOptions[2] = activeDobNode.val() + " " + activeTobNode.val();
        GlobalCurrentSaveOptions[3] = activeAccNode.is(':checked') | 0;
        /* PrintArray(GlobalCurrentSaveOptions);  */
        _f__populateHMSArraysFromValueArraysValues();
        _f__populateValueArraysFromScreenValues();
        _f__InitProcessedDataArrayToBlanksIfUndefined();
        if (IsSubjectSave) {
            if (!_f__RecordEntriesAreNoneZero(SubjectsValues) || activeNameNode.val().length < 1) {
                _f__warnClear();
                _f__warnP("Enter some values, Choose your subjects name");
                activeNameNode.val("Subject");
                return;
            }
            _f__LoadSubjectChartData(activeNameNode.val());
        } else {
            if (!_f__RecordEntriesAreNoneZero(TargetsValues) || activeNameNode.val().length < 1) {
                _f__warnClear();
                _f__warnP("Enter some values, Choose your " + NameOfColumn + " name");
                activeNameNode.val("Target");
                return;
            }
            _f__LoadTargetChartData(activeNameNode.val());
        }
        _f__populateHMSArraysFromValueArraysValues();
        _f__populateValueArraysFromScreenValues();
        /* 	//debugP("GlobalCurrentSaveOptions in the save path");  */
        _f__doConditionalSave(NameOfColumn, activeNameNode.val(), GlobalCurrentSaveOptions);
    }
function _f__doConditionalSave(ColumnToSave, nameOfRecord, SaveOptions) {
        var positionsToSave; // array
        var arrayToSave; // Theme
        if (ColumnToSave == "subject") {
            arrayToSave = SubjectProcessedData;
            positionsToSave = SubjectsValues;
        } else {
            if (ColumnToSave == "target") {
                arrayToSave = TargetsProcessedData;
                positionsToSave = TargetsValues;
            }
        }
        /*//debugP(positionsToSave);*/
        if (!localStorage.getItem(nameOfRecord)) {
            _f__StorageSaveOperation(nameOfRecord, positionsToSave, arrayToSave, SaveOptions);
        } else {
            _f__deleteRecordFromRecordsList(nameOfRecord);
            _f__StorageSaveOperation(nameOfRecord, positionsToSave, arrayToSave, SaveOptions);
        }
    }

function setSelectedIndex(s, v) {

    for ( var i = 0; i < s.options.length; i++ ) {

        if ( s.options[i].text == v ) {

            s.options[i].selected = true;

            return;

        }

    }

}

function _f__ClickedLoadForImport(ColumnToSave) {
        var placeToLoad;
        if (ColumnToSave.toLowerCase() == "subject") {
            placeToLoad = $("#SubjectName").prop('value'); // SubjectName
        } else {
            placeToLoad = $("#TargetName").prop('value'); // SubjectName		
        }


        recordsList = document.getElementById('person');

		setSelectedIndex(recordsList,placeToLoad);

        
		var index = $('#person option:selected');
        if (index.length > 1) {
            return;
        }
        if (index.length == 1) {
            _f__populateValueArraysFromStoredValues($("#person option:selected").text(), ColumnToSave);
        } else {
        }
        _f__populateHMSArraysFromValueArraysValues();
    }


function _f__ClickedLoadForTests(ColumnToSave) {
        var placeToLoad;
        if (ColumnToSave.toLowerCase() == "subject") {
            placeToLoad = $("#SubjectName").prop('value'); // SubjectName
        } else {
            placeToLoad = $("#TargetName").prop('value'); // SubjectName		
        }
        recordsList = document.getElementById('person');
        var index = $('#person option:selected');
        if (index.length > 1) {
            return;
        }
        if (index.length == 1) {
            _f__populateValueArraysFromStoredValues($("#person option:selected").text(), ColumnToSave);
        } else {
        }
        _f__populateHMSArraysFromValueArraysValues();
    }
function _f__ClickedLoad(event, ColumnToSave) {
        var placeToLoad;
        if (ColumnToSave == "subject") {
            placeToLoad = $("#SubjectName").prop('value'); // SubjectName
        } else {
            placeToLoad = $("#TargetName").prop('value'); // SubjectName		
        }
        recordsList = document.getElementById('person');
        var index = $('#person option:selected');
        if (index.length > 1) {
            _f__warnClear();
            _f__warnP("Please select a single record from the list to load");
            return;
        }
        if (index.length == 1) {
            _f__populateValueArraysFromStoredValues($("#person option:selected").text(), ColumnToSave);
        } else {
            _f__warnClear();
            _f__warnP("Please select a single record from the list to load");
        }
        _f__populateHMSArraysFromValueArraysValues();
        _f__clearCharts();
    }
function _f__populateValueArraysFromScreenValues() {
        var colNum = 12;
        var HorX = new Array();
        var j = 0;
        var i = 0;
        for (i = 0; i < colNum; i++) {
            if (document.f1.numvalue1[i].value != '' && parseInt(document.f1.numvalue1[i].value) >= 0) {
                var PlanetNeatTitles = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Ascendant", "Midheaven"];
                var PlanetNeatSymbols = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"];
                HorX[j] = PlanetNeatSymbols[i]; //   document.f1.label[i].value;
                SubjectsValues[j] = document.f1.numvalue1[i].value;
                TargetsValues[j] = document.f1.numvalue2[i].value;
                j++;
            }
        }
    }
function _f__CreateOptionsObject() {
        _f__CreateOptionsObjectWithAssumedState(); // need this here for xProfile
        _f__SetOptionsObjectFromWidgetState();
        /* 	PrintArray(OptionsObject); */
        _f__PopulatePerUserRecordVariableFromPassedPerlVariablesOrScreenValues();
        /* 	PrintArray(GlobalCurrentSaveOptions); */
        /* 	PrintArray(OptionsObject);  */
        _f__AppendResearchFlagsToGlobalOptionsObject();
        /*  //debugP("final GlobalCurrentSaveOptions after _f__CreateOptionsObject"); */
        /* PrintArray(GlobalCurrentSaveOptions);  */
    }
function _f__CreateOptionsObjectButPreservePerUserRecordVariablesToCurrentState() {
    _f__SetOptionsObjectFromWidgetState();
    /* 	PrintArray(GlobalCurrentSaveOptions); */
    _f__AppendResearchFlagsToGlobalOptionsObject();
    /* 	PrintArray(GlobalCurrentSaveOptions);  */
}

function _f__AppendPreparedStateNotSetFromWidgetsThenCreateOptions(NameOfColumn) {
        _f__PopulatePerUserRecordVariableFromPassedPerlVariablesOrScreenValues(NameOfColumn);
        _f__AppendResearchFlagsToGlobalOptionsObject();
    }
function _f__PopulatePerUserRecordVariableFromPassedPerlVariablesOrScreenValues(NameOfColumn) {
        var IsSubjectSave = false;
        if (typeof NameOfColumn === "undefined") {
            IsSubjectSave = true;
            /* alert ("creating a default object..."); */
        } else if (NameOfColumn == "subject") {
            IsSubjectSave = true;
        }
        if (parseFloat(lat).IsDefault()) {
            if (IsSubjectSave) {
                lat = $('#uk-lat').val().trim();
            } else {
                lat = $('#uk-lat-tar').val().trim();
            }
        };
        if (parseFloat(lon).IsDefault()) {
            if (IsSubjectSave) {
                lon = $('#uk-lon').val().trim();
            } else {
                lon = $('#uk-lon-tar').val().trim();
            }
        };
        var activeDateNode = $('#uk-dob'); // Defaults to subject!
        if (dob.IsDefaultDOB()) {
            if (IsSubjectSave) {
                activeDateNode = $('#uk-dob');
            } else {
                activeDateNode = $('#uk-dob-tar');
            }
            /* _f__warnP("using undefined dob path - since this is all blanks"); */
            dob = activeDateNode.val().trim().replace(/\//g, ' ');
            if (!$.isNumeric(dob.replace(/ /g, ''))) {
                _f__warnP(dob + " doesnt look numeric...")
            }
        };
        if (dob.trim() == __DOB_DEFAULTS__) {
            if (activeDateNode.val().trim().replace(/\//g, ' ') == __DOB_DEFAULTS__) {
            } else {
                dob = activeDateNode.val().trim().replace(/\//g, ' ');
                if (!$.isNumeric(dob.replace(/ /g, ''))) {
                    _f__warnP(dob + " doesnt look numeric...");
                }
            }
        }
        else {
            if (dob.trim() != activeDateNode.val().trim().replace(/\//g, ' ')) {
                if (activeDateNode.is(':visible')) {
                    if (activeDateNode.val().trim().replace(/\//g, ' ') != __DOB_DEFAULTS__) {
                        dob = activeDateNode.val().trim().replace(/\//g, ' ');
                    }
                }
            }
        }
        var activeTimeNode = $('#uk-tob');
        if (tob.IsDefaultTOB()) {
            if (IsSubjectSave) {
                activeTimeNode = $('#uk-tob');
            } else {
                activeTimeNode = $('#uk-tob-tar');
            }
            tob = activeTimeNode.val().trim().replace(/:/g, ' ');
            if (!$.isNumeric(tob.replace(/ /g, ''))) {
                _f__warnP(tob + " doesnt look numeric...")
            }
        };
        if (tob.trim() == __DOB_DEFAULTS__) {
            if (activeTimeNode.val().trim().replace(/\//g, ' ') == __TOB_DEFAULTS__) {
            } else {
                tob = activeTimeNode.val().trim().replace(/\//g, ' ');
                if (!$.isNumeric(tob.replace(/ /g, ''))) {
                    _f__warnP(tob + " doesnt look numeric...");
                }
            }
        }
        else {
            if (tob.trim() != activeTimeNode.val().trim().replace(/\//g, ' ')) {
                if (activeTimeNode.is(':visible')) {
                    if (activeTimeNode.val().trim().replace(/\//g, ' ') != __TOB_DEFAULTS__) {
                        tob = activeTimeNode.val().trim().replace(/\//g, ' ');
                    }
                }
            }
        }
        var activeAccNode = $('#uk-acc');
        if (parseInt(acc).IsDefault()) {
            if (IsSubjectSave) {
                activeAccNode = $('#uk-acc')
            } else {
                activeAccNode = $('#uk-acc-tar');
            }
            acc = 1 * (activeAccNode.is(':checked'));
        } else
        if (activeAccNode.is(':visible')) {
            if (acc != 1 * (activeAccNode.is(':checked'))) {
                acc = 1 * (activeAccNode.is(':checked'));
            }
        }
        if (lat.length == 0) {
            lat = "0";
        };
        if (lon.length == 0) {
            lon = "0";
        };
    }
function _f__CreateOptionsObjectWithAssumedState() {
    OptionsObjectLon = lon;
    OptionsObjectLat = lat;
    OptionsObjectDateTimeAccurate = acc;
    OptionsObjectDateTime = "";
    GlobalCurrentSaveOptions = [
        OptionsObjectLon, // LONG
        OptionsObjectLat, // LAT
        OptionsObjectDateTime, // DATETIME of BIRTH
        OptionsObjectDateTimeAccurate, // Certain of TOB
        OptionsObjectUtOffset, //OptionsObjectUtOffset, // seconds
        OptionsObjectDstOffset //OptionsObjectDstOffset, // seconds
    ];
}

function _f__AppendResearchFlagsToGlobalOptionsObject() {
        $.each(OptionsObject, function (index, value) {
            GlobalCurrentSaveOptions[GlobalCurrentSaveOptions.length] = value;
        });
        /*  $.each(GlobalCurrentSaveOptions, function  (index, value) {
	 

	 })
	 */
    }
function _f__buildOptions() {
        var listOfIdsToEnable = [];
        var htmlAppend = "";
        jQuery.each(OptionObjectValues, function (key, value) {
            var outer = key;
            var arry = value;
            htmlAppend += '<span class="span-radio"><span class="span-radio-title">' + outer.toLowerCase().substring(0, 10) + '</span><br/>';
            if (arry[0] == -1) {
                debugP("YAY");
                var command = 'javascript:OptionsObject.' + outer + '=' + this.text;
                htmlAppend += '	  <input class="options-radio" style="height: 20px;" onchange="' + command + ' ;_f__DoCompleteOptionsSaveThenLoadCycle();" type="text" id=' + outer + '' + value + ' name="' + outer + '" value="' + arry[1] + '"/><br/>';
                listOfIdsToEnable.push([(outer + '' + value), 'Set' + '' + outer + '(' + value + ')']);
            } else {
                jQuery.each(arry, function (key, value) {
                    var command = 'javascript:OptionsObject.' + outer + '=' + value;
                    htmlAppend += '	  <input class="options-radio" onchange="' + command + ' ;_f__DoCompleteOptionsSaveThenLoadCycle();" type="radio" id=' + outer + '' + value + ' name="' + outer + '" value="' + value + '"/>' + key + '<br/>';
                    listOfIdsToEnable.push([(outer + '' + value), 'Set' + '' + outer + '(' + value + ')']);
                })
            }
            htmlAppend += "</span>";
        })
        $('#optionsHolder').html(htmlAppend);
        _f__GetGlobalOptions(OptionsObject);
    }
function _f__InitProcessedDataArrayToBlanksIfUndefined() {
    if (typeof SubjectProcessedData === "undefined") {
        SubjectProcessedData = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    }
    if (typeof TargetsProcessedData === "undefined") {
        TargetsProcessedData = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    }
}

function _f__LoadBothChartsData(subjectName, targetName) {
        _f__LoadSubjectChartData(subjectName);
        _f__LoadTargetChartData(targetName);
    }
function _f__LoadSubjectChartData(name) {
        _f__SetNativityYear(name);
        getThemeValues(SubjectsValues[0], SubjectsValues[1], SubjectsValues[2], SubjectsValues[3], SubjectsValues[4], SubjectsValues[5], SubjectsValues[6], SubjectsValues[7], SubjectsValues[8], SubjectsValues[9], SubjectsValues[10], SubjectsValues[11]);
        SubjectProcessedData = theme.slice();
    }
function _f__LoadTargetChartData(name) {
    _f__SetNativityYear(name);
    getThemeValues(TargetsValues[0], TargetsValues[1], TargetsValues[2], TargetsValues[3], TargetsValues[4], TargetsValues[5], TargetsValues[6], TargetsValues[7], TargetsValues[8], TargetsValues[9], TargetsValues[10], TargetsValues[11]);
    TargetsProcessedData = theme.slice();
}

function _f__UpdateAcc() {
        acc = 1 * ($('#uk-acc').is(':checked'));
    }
function _f__storageAvailable(type) {
        try {
            var storage = window[type],
                x = '__storage_test__';
            storage.setItem(x, x);
            storage.removeItem(x);
            return true;
        } catch (e) {
            return false;
        }
    }
function _f__warnDefaulted(name, value) {
        _f__warnP("defaulted " + name + " to " + value);
    }
function _f__ToggleOptionsVisiblity() {
        $('. PreviouslyHideableOptionsComponent').toggleClass("options-hidden-record");
    }
function _f__DebugAlertOptions() {
        /*	//debugP("OPTIONS: aoIndex:"+OptionsObject.aoIndex+"tfIndex:"+OptionsObject.tfIndex+"orbType:"+OptionsObject.orbType+"precessionFlag:"+OptionsObject.precessionFlag   );	*/
    }
function _f__SaveChanges() {
        /* 	$('#uk-dob').val(dayOB+"/"+monOB+"/"+yearOB);
        	$('#uk-tob').val() */
        $('#uk-dialog').modal('hide');
        warnP("TODO: Make this save just the DOB/TOB  LONG/LAT when offline");
    }
function _f__Cleanse(thing) {
    var ret = "";
    var forbiddenChars = new RegExp("[^a-zA-Z0-9]", 'g');
    if (forbiddenChars.test(thing)) {
        ret = thing.replace(forbiddenChars, '');
    }
    return ret;
}

function PrintArray(thing) {
    var stringToPrint = "";
    $.each(thing, function (index, value) {
        /* stringToPrint += index +" is "+ value+"\n"; */
        stringToPrint += index + " is " + value + "<br />";
    })
    _f__infoP(stringToPrint);
}
Object.size = function (obj) {
    var size = 0,
        key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};

function _f__buildStringForInfo(theString) {}

function _f__AppendRecordHTML() {
    var externalIndexRef = 0;
    var internalModuloIndex = 1;
    var modulus = 12;
    $(".recordHolder").each(function () {
        /*      $(this).html("<tr class='record'><td class='record'><span style='white-space: nowrap'><span class='ui-icon  ui-icon-search sign-selector' style='display:inline-block; '></span><select  name='SPI' class='SignPrefix' ><option class='SignPrefixChoice' value='0'>Aries</option><option class='SignPrefixChoice' value='30'>Taurus</option><option class='SignPrefixChoice' value='60'>Gemini</option><option class='SignPrefixChoice' value='90'>Cancer</option><option class='SignPrefixChoice' value='120'>Leo</option><option class='SignPrefixChoice' value='150'>Virgo</option><option class='SignPrefixChoice' value='180'>Libra</option><option class='SignPrefixChoice' value='210'>Scorpio</option><option class='SignPrefixChoice' value='240'>Sagittarius</option><option class='SignPrefixChoice' value='270'>Capricorn</option><option class='SignPrefixChoice' value='300'>Aquarius</option><option class='SignPrefixChoice' value='330'>Pisces</option></select><input class='hmsBox' type='text' size='20'  name='hms' value='' >&nbsp;<input class='numVals hideable' type='decimal' size='4'  name='numvalue" + internalModuloIndex + "' value='0' disable onblur='ConvertDecimalToHMSSingle(this)'><div class='sign-label unset-sign'>ARI</div></span></td></tr>");
         */
        var prettyColName = "";
        if (internalModuloIndex == 1) {
            prettyColName = "subject"
        } else {
            prettyColName = "target"
        }
        $(this).html("<tr class='record'><td class='record'><span style='white-space: nowrap'><span class='ui-icon ui-icon-search sign-selector' style='display:inline-block; '></span><select id='combobox' class='SignPrefix'> <option class='SignPrefixChoice' value='0'>Aries</option><option class='SignPrefixChoice' value='30'>Taurus</option><option class='SignPrefixChoice' value='60'>Gemini</option><option class='SignPrefixChoice' value='90'>Cancer</option><option class='SignPrefixChoice' value='120'>Leo</option><option class='SignPrefixChoice' value='150'>Virgo</option><option class='SignPrefixChoice' value='180'>Libra</option><option class='SignPrefixChoice' value='210'>Scorpio</option><option class='SignPrefixChoice' value='240'>Sagittarius</option><option class='SignPrefixChoice' value='270'>Capricorn</option><option class='SignPrefixChoice' value='300'>Aquarius</option><option class='SignPrefixChoice' value='330'>Pisces</option></select><input class='hmsBox hms-" + prettyColName + "' type='text' size='20'  name='hms' value='' >&nbsp;<input class='numVals numVals-" + prettyColName + " hideable' type='decimal' size='4'  name='numvalue" + internalModuloIndex + "' value='0' disable onblur='_f__ConvertDecimalToHMSSingle(this)'><div class='sign-label unset-sign'>ARI</div></span></td></tr>");
        externalIndexRef++;
        if (externalIndexRef >= modulus) {
            internalModuloIndex++;
            externalIndexRef = 0;
        }
    });
    $(".seconds, .minutes, .hours, .numVals, .hmsBox").each(function () {
        $(this).uniqueId()
    });
}

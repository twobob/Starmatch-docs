var DOING_MAP_LOOKUP=false;

var REPLACE_LETTERS_WITH_SPACE_LENGTH = 10;


var time_delay =  2500;


var INPUT_BOX_LENGTH = 17 // old value 8;
var INPUT_LAT_LON_LENGTH = 17; // old value 6;

var CURRENT_LOCATION_LAT_LON = 17   // old value 8;
var LAT_LON_MAXLENGTH = 17; // old value 6;


var FILE_LOAD_TIMEOUT = 400;
var IMPORT_POPUP_TIMEOUT = 10;

var SubjectProcessedData;
var TargetsProcessedData;
var SubjectsValues = new Array();
var TargetsValues = new Array();
var response = "";

var oldUkLatValue = 0;
var oldUkLonValue = 0;

const __LON_DEFAULTS__ = 0;
const __LAT_DEFAULTS__ = 0;


const __DOB_DEFAULTS__ = "00 00 0000";
const __TOB_DEFAULTS__ = "00 00";

const __NUMBER_DEFAULTS__ = Number.MIN_VALUE;
const __STRING_DEFAULTS__ = "__INIT_ME__";


var dob = __DOB_DEFAULTS__;
var tob = __TOB_DEFAULTS__;
var lon = __NUMBER_DEFAULTS__;
var lat = __NUMBER_DEFAULTS__;
var acc = __NUMBER_DEFAULTS__;

var command = __STRING_DEFAULTS__;
var datastream = __STRING_DEFAULTS__;
var datacommand = __STRING_DEFAULTS__;
var dataPassed = __STRING_DEFAULTS__;

String.prototype.IsDefaultTOB = function () {
    return this == __TOB_DEFAULTS__;
};

String.prototype.IsDefaultDOB = function () {
    return this == __DOB_DEFAULTS__;
};

String.prototype.IsDefault = function () {
    return this == __STRING_DEFAULTS__;
};

Number.prototype.IsDefault = function () {
    return this == __NUMBER_DEFAULTS__;
};

var CRLF = "&#013;&#010;";
var br = "<br />";
var originalentries = null;
var recordsList; // = document.getElementById('person');

var latsign = 1;
var absdlat = 0;
var absmlat = 0;
var absslat = 0;

var SignNeatTitles = ["Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo", "Libra", "Scorpio", "Sagittarius", "Capricon", "Aquarius", "Pisces"];

var useLookups = 1;
var maphtml = "";

/*
function _f__InitTheFloatingVariables() {
}
*/
function _f__SetGlobalOptions(reference) {
    _f__CreateOptionsObject(reference);
    _f__SaveOptionsCookie(reference);
}

function _f__SaveOptionsCookie(reference) {
    Cookies.set('OptionsObject', reference, {
        expires: 8364
    });
}

function _f__DoCompleteOptionsSaveThenLoadCycle() {
    _f__SetGlobalOptions(OptionsObject);
    _f__GetGlobalOptions(OptionsObject);
}

function _f__GetGlobalOptions(reference) {
    reference = Cookies.getJSON('OptionsObject');

    if (typeof reference === "undefined") {

        reference = OptionsObject; // defaults
        _f__SetGlobalOptions(reference); // save it back
    }

    /* PrintArray(reference); */

    _f__SetWidgetStateFromOptionsObject(reference);

}

var cc = {
    cc: "false"
};
var ccToast = "";

function _f__CheckCcCookieOptions(reference) {
    _f__GetCcCookie(cc);



    if (cc['cc'] == "true")
        return;


    toastr.clear();

    toastr.options = {
        "closeButton": false,
        "debug": false,
        "newestOnTop": false,
        "progressBar": false,
        "positionClass": "toast-top-right-cc",
        "preventDuplicates": false,
        "showDuration": "300",
        "hideDuration": "1000",
        "timeOut": 0,
        "extendedTimeOut": 0,
        "showEasing": "swing",
        "hideEasing": "linear",
        "showMethod": "fadeIn",
        "hideMethod": "fadeOut",
        "tapToDismiss": false
    }

    toastr.options.onclick = _f__CookieCcCallback;

    ccToast = toastr["cc"]("We may place cookies on your computer to store your options. <br/>More details can be found in our privacy policy &nbsp;<button type='button' class='clear cc-ok btn clear'>Okay, I get it.</button>", "COMPLIANCE");

}

function _f__CookieCcCallback() {
    ccToast.remove();



    cc = {
        cc: "true"
    };

    _f__SaveCcCookie(cc);
    _f__CheckCcCookieOptions(cc);

}

function _f__GetCcCookie(reference) {
    reference = Cookies.getJSON('cc');

    if (typeof reference === "undefined") {
        reference = {
            cc: "false"
        }; // defaults
        _f__SaveCcCookie(reference); // save it back
    }

    cc = reference;
}

function _f__UndoCc() {

    cc = {
        cc: "false"
    };
    _f__SaveCcCookie(cc);

}

function _f__SaveCcCookie(reference) {
    Cookies.set('cc', reference, {
        expires: 2147
    });
}

function _f__DMS2Decimaled(days, minutes, seconds) {
	 var decimaled = "Error?" 
	 if(days < 0) 
	 {    decimaled = parseFloat(days) -     (  parseFloat(minutes)/60) - (parseFloat(seconds)/(60*60));	 }
	 else
	 {    decimaled = parseFloat(days) +     (  parseFloat(minutes)/60) + (parseFloat(seconds)/(60*60));	 } 	
	return decimaled;
}

function _f__LimitCertainKeystrokes() {
    $('#uk-lat,#uk-lat-tar').focus(function () {
        oldUkLatValue = $(this).val();
        $(this).val('');
    }).blur(function () {
        if ($(this).val() == '') {
            $(this).val(oldUkLatValue);
        } else {
            $(this).val($(this).val().replace(/[^\d.-]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, '').trim());
        }
		
		var dms = $(this).val().split(" ");
		var spaceCount = (dms.length - 1)
		if (spaceCount ==2) {$(this).val(_f__DMS2Decimaled(dms[0],dms[1],dms[2]  )	)  }
		if (spaceCount ==1) {$(this).val(_f__DMS2Decimaled(dms[0],dms[1],0  )	)  }
		
        lat = $(this).val();q
    })[0].maxLength = INPUT_BOX_LENGTH;

    $('#uk-lon,#uk-lon-tar').focus(function () {
        oldUkLonValue = $(this).val();
        $(this).val('');
    }).blur(function () {
        if ($(this).val() == '') {
            $(this).val(oldUkLonValue);
        } else
			{
            $(this).val($(this).val().replace(/[^\d.-]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, ''));
        }
		var dms = $(this).val().split(" ");
		var spaceCount = (dms.length - 1)
		if (spaceCount ==2) {$(this).val(_f__DMS2Decimaled(dms[0],dms[1],dms[2]  )	)  }
		if (spaceCount ==1) {$(this).val(dms[0]) ; _f__infoP("DECIMAL or D M S" ,2000); }
        lon = $(this).val();
    })[0].maxLength = INPUT_BOX_LENGTH;

    $('#uk-picker').locationpicker({
        location: {
            latitude: 51.566520,
            longitude: -0.1394648
        },
        radius: 0,
        inputBinding: {
            latitudeInput: $('#uk-lat'),
            longitudeInput: $('#uk-lon'),
            locationNameInput: $('#uk-address')
        },
        enableAutocomplete: true,
        onchanged: function (currentLocation) {
        }
    });
    $('#uk-dialog').on('shown.bs.modal', function () {
        $('#uk-picker').locationpicker('autosize');
    });



    $('#uk-dob, #uk-dob-tar').bind('focus', function () {
            if ($(this).val() == '00 00 0000') {
                $(this).val('');
            }
        })
        .bind('blur', function () {
            if ($(this).val().trim() == '') {
                $(this).val('00 00 0000');
            }
            $(this).val($(this).val().replace(/[^\d.-]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, ''));

            var ukTobVal = $('#uk-tob').val();

            var splitDob = $(this).val().split(' ');
            var splitTob = ukTobVal.split(' ');


            if ($(this).val().trim() != '' && ukTobVal != '00 00' && ukTobVal != '') {
                OptionsObjectDateTime = new Date(splitDob[2], splitDob[1], splitDob[0], splitTob[0], splitTob[1]);
            }


        })


        .on('paste', function () {

            var self = $(this);
            setTimeout(function () { // Magically replace letters with a space.
                $(self).val($(self).val().replace(/[^\d.-]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, ''));
            }, REPLACE_LETTERS_WITH_SPACE_LENGTH);

        }).on('keyup', function () {

            var self = $(this);
            setTimeout(function () { // Magically replace letters with a space.
                $(self).val($(self).val().replace(/[^\d.-]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, ''));
            }, REPLACE_LETTERS_WITH_SPACE_LENGTH);

        })

    ;
    $('#uk-tob, #uk-tob-tar').bind('focus', function () {
            if ($(this).val() == '00 00') {
                $(this).val('');
            }
        })
        .bind('blur', function () {
            if ($(this).val().trim() == '') {
                $(this).val('00 00');
            }
            $(this).val($(this).val().replace(/[^\d.]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, ''));


            var ukDobVal = $('#uk-dob').val();

            var splitTob = $(this).val().split(' ');
            var splitDob = ukDobVal.split(' ');


            if ($(this).val().trim() != '' && ukDobVal != '00 00 00' && ukDobVal != '') {
                OptionsObjectDateTime = new Date(splitDob[2], splitDob[1], splitDob[0], splitTob[0], splitTob[1]);
            }


            if ($(this).val().trim() != '') {
                OptionsObjectDateTime = $(this).val();
            }
        })
        .on('paste', function () {
            var self = $(this);
            setTimeout(function () { // Magically replace letters with a space.
                $(self).val($(self).val().replace(/[^\d.]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, ''));
            }, REPLACE_LETTERS_WITH_SPACE_LENGTH);

        }).on('keyup', function () {

            var self = $(this);
            setTimeout(function () { // Magically replace letters with a space.
                $(self).val($(self).val().replace(/[^\d.]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, ''));
            }, REPLACE_LETTERS_WITH_SPACE_LENGTH);

        });

    /*} */

}


function initToZero(nameOfThing){
	
	$('#'+nameOfThing).val('0');
	
}

function _f__InitPage() {
    useLookups = 1;
	
	initToZero('uk-lat-tar');
	initToZero('uk-lon-tar');
	initToZero('uk-lon');
	initToZero('uk-lat');
	
	
    new Clipboard('.copy-btn');

    var fileInput = document.getElementById('fileInput');
	
	 var coeffsInput = document.getElementById('CoeffsInput');
	
    var fileDisplayArea = document.getElementById('fileDisplayArea');
	
	
	function handleFileSelect(event) {
    if (window.File && window.FileList && window.FileReader) {

        var files = event.target.files; //FileList object
        var output = document.getElementById("result");

		
		if(files.length != 2)
		{
			_f__infoP("REMEMBER TO SELECT TWO FILES", 1000);
		}
		
		   var ourNamesFileReader = new FileReader();
			
			var ourDataFileReader = new FileReader();

		   var textType = /text.*/;
		
        for (var i = 0; i < files.length; i++) {
            var file = files[i];

			
			if(file.name.includes("_names.txt"))
			{
			
            ourNamesFileReader.onload = (function(e) 
            {
				
				 var contents = ourNamesFileReader.result;
               curveNameSet =    JSON.parse(contents).split("\n");
            });

            ourNamesFileReader.readAsText(file);
			
			
			 }
			
			if(file.name.includes("_data.txt"))
			{
			
            ourDataFileReader.onload = (function(eData) 
            {
				 var contents = ourDataFileReader.result;
              curveDataSet = JSON.parse(contents);
            });

            ourDataFileReader.readAsText(file);
			 }
		}
    } else {
        console.log("Your browser does not support File API");
    }
	setTimeout(_f__coeffAnalysis, FILE_LOAD_TIMEOUT);
	setTimeout(	CloseTheImportWindow, IMPORT_POPUP_TIMEOUT);
}

function CloseTheImportWindow(){
		$( "#exporterHolder" ).dialog('close');
}



coeffsInput.addEventListener('change', handleFileSelect, false);
	
	

var totalDelay =0;

function setup_reader(files, i) {



if (typeof window.FileReader !== 'function') {
        alert("The file API isn't supported on this browser yet.");
    }

    var file = files[i];
    var name = file.name;
    var reader = new FileReader();
    reader.onload = function(e){
                        readerLoaded(e, files, i, name);
                    };
    reader.readAsBinaryString(file);
}



function readerLoaded(e, files, i, name) {
    


    var bin = e.target.result;

		setTimeout(function(){ 
   
	debugP('importing '+files[i].name);
   	_f__handle_returnIMPORT(bin, files[i].name);

}, totalDelay);

	totalDelay += (bin.length / 100.0) * time_delay;
    if (i < files.length - 1) {
        setup_reader(files, i+1);
    }else{


 debugP('processing delay '+(totalDelay/ 1000).toPrecision(5) );


setTimeout(function(){ 
   
	_f__infoP("IMPORTS COMPLETE", 5000);

}, totalDelay+ 5000);
}
}



   fileInput.addEventListener('change', function (e) {
totalDelay = 0;
  setup_reader(fileInput.files, 0);



});

	var nameSelectionListIMPORT = [];


function _f__handle_returnIMPORT(returnedString, filenameBeingProcessed) {
	

	nameSelectionListIMPORT = [];
	ClearCurrentSelection();
	
	
	var updateCount = 0;
	
	
    var uncoded = returnedString.trim(-1).split(delimeter);

    var name = "";

	var nameList = "";
	
	
	
    $.each(uncoded, function (index, val) {

        var result = "";

		if (val == "")
		{return;}
		
		
	
	
	 name = val.split('|')[0];
	
	
	var dataGetter = val.split('|')[1];
	
        var secondLeveluncoded = dataGetter.split(',');


		
		updateCount ++;
		
	
	var myNewArray = dataGetter.split(',');
		
		
		var SubjectsValues =  [myNewArray[0], myNewArray[2], myNewArray[4],myNewArray[6], myNewArray[8],myNewArray[10], myNewArray[12], myNewArray[14], myNewArray[16],myNewArray[18], myNewArray[20], myNewArray[22]     ];
		
		var subjectProcessedData =  [myNewArray[1], myNewArray[3], myNewArray[5],myNewArray[7], myNewArray[9],myNewArray[11], myNewArray[13], myNewArray[15], myNewArray[17],myNewArray[19], myNewArray[21], myNewArray[23]     ];
		
		
		
		var myShakyOptionsArray = myNewArray.slice(24);
		
		name = name.replace('--',', ');
		name = name.replace(',-',', ');
		
		
		
		
		nameSelectionListIMPORT[nameSelectionListIMPORT.length] =  name;
		
		
		var addon = ", ";
		if(updateCount % 3 == 0)
		{
			addon = " ...<br />";
		}
		nameList += "["+updateCount+": "+ name+"] "+addon;
		
		_f__StorageSaveOperation(name,SubjectsValues,subjectProcessedData,myShakyOptionsArray   );

    })
	
	_f__updatedP("Updated "+updateCount+" Records for<br />"+nameList);
	
	
	


	 $.each(nameSelectionListIMPORT, function (index, val){
		 
setTimeout(function(){
	
		$('#SubjectName').val(val);
		
		_f__ClickedLoadForImport('subject');

	
		 _f__SaveChangesGetChart('subject', false, val);

		}, (index * time_delay) +2);

		 
	 });

	
	
	setTimeout(function(){__f_SelectAndSaveIMPORT(filenameBeingProcessed) }, (nameSelectionListIMPORT.length * (time_delay + 200))  );
	
	
	
}


function __f_SelectAndSaveIMPORT(filenameBeingProcessed){
	__f_SelectListIMPORT(nameSelectionListIMPORT); 
	SaveCurrentSelectionWithName(filenameBeingProcessed.substr(0, filenameBeingProcessed.length-4));
	
	SaveCurrentSelectionWithName(filenameBeingProcessed.substr(0, filenameBeingProcessed.length-4)+' TEST')
	
}


function __f_SelectListIMPORT(listToSelect){
	
	$(listToSelect).each( function(ind, val){     __f_SelectbyNameIMPORT(val)    }  )
	
}


	

function  __f_SelectbyNameIMPORT(nameToSelect)
{
	
	$('#person option').each( function(ind, val){  if (  $(val).text() == nameToSelect   ){    $(val).prop('selected', true)     }                        }  )
	
}

	
	
	
	
	
	

	

    _f__GetGlobalOptions(OptionsObject);

    _f__buildOptions();

    $("#optionsHolder").dialog({
        autoOpen: false,
        width: 180,
        minWidth: 180,
        height: 600,
        minHeight: 600,
        modal: false,
        resizable: true,
        title: "Options",
        dialogClass: "optionsHolder",
        position: {
            my: "left bottom",
            at: "right top",
            of: "#logo"
        },
        show: {
            effect: "blind",
            duration: 300
        },
        hide: {
            effect: "blind",
            duration: 300
        }

    });


    $("#chartsHolder").dialog({
        autoOpen: false,
        width: 400,
        minWidth: 400,
        height: 720,
        minHeight: 720,
        modal: false,
        resizable: true,
        title: "Charts",
        dialogClass: "chartsHolder",
        position: {
            my: "left bottom",
            at: "right top",
            of: "#logo"
        },
        show: {
            effect: "blind",
            duration: 300
        },
        hide: {
            effect: "blind",
            duration: 300
        },
        buttons: {
            "Report": function () {

                lastFullChartReport.sort(function (a, b) {
                    return b[1] - a[1];
                });


                _f__flagP(lastFullChartReport);
            },
            "Close": function () {
                $(this).dialog("close");
            }
        }

    });


    $("#xchartsHolder").dialog({
        autoOpen: false,
        width: 400,
        minWidth: 400,
        height: 630,
        minHeight: 610,
        modal: false,
        resizable: true,
        title: "xProfile Charts",
        dialogClass: "xchartsHolder",
        position: {
            my: "left bottom",
            at: "right top",
            of: "#logo"
        },
        show: {
            effect: "blind",
            duration: 300
        },
        hide: {
            effect: "blind",
            duration: 300
        },
        buttons: [{
                text: "Report",
                click: function () {
                    _f__flagP(lastFullChartReport);
                },
                class: "overideReport",
                style: "color:#fff !important"
            },
            {
                text: "Close",
                click: function () {
                    $(this).dialog("close");
                },
                class: "overideClose",
                style: "color:#fff !important"
            }
        ]

    });



    $("#exporterHolder").dialog({
        autoOpen: false,
        width: 400,
        minWidth: 400,
        height: 280,
        minHeight: 200,
        modal: false,
        resizable: false,
        title: "EXPORT / IMPORT",
        dialogClass: "exporterHolder",
        position: {
            my: "left bottom",
            at: "right top",
            of: "#logo"
        },
        show: {
            effect: "blind",
            duration: 300
        },
        hide: {
            effect: "blind",
            duration: 300
        }
        /*   ,
	  buttons: {
        "Save Current Selection": function() {
			
			
        SaveCurrentSelection();
        },
        "Close": function() {
          $( this ).dialog( "close" );
        }
	  } */

    });

    _f__BindKeys();
    _f__AppendRecordHTML();
    _f__populateRecordsList();


    _f__populateSelectionsList();


    _f__InitSearching();

    _f__populateValueArraysFromScreenValues();

    _f__populateHMSArraysFromValueArraysValues();
	
   _f__LimitCertainKeystrokes();

    _f__AutoSetupTheTabbingOrder();

    /*        $('#clearlasttoast').click(function  () {
            toastr.clear(getLastToast());
        }); */

    $('#cleartoasts').click(function () {

        toastr.clear();
        $(this).blur();
    });


    if (typeof (dataPassed) === "undefined") {
        useLookups = 0;
    } else {
        var maphtml = "<div class='modal-dialog' style='background-color: rgb(202, 222, 230);'> <div class='modal-content' style='background-color: rgb(202, 222, 230);'> <div class='modal-header' style='background-color: rgb(202, 222, 230);'> <button type='button' class='close' data-dismiss='modal' aria-label='Close'><span aria-hidden='true'>&times;</span></button> <h4 class='modal-title'>Enter Your Birth Details</h4> </div> <div class='detailsForm' style='pointer-events: none; height: 130px; margin-top: 10px;' > <form id='f2' style='width:800px;' name='f2' class='n-t-small' action='index.pl' method='get'> <div class='entry-holder' > <span style='white-space:nowrap;'> <label style='padding-right: 20px; width:180px ' for='uk-dob-map' class='control-label defaultedValueDob' >DOB</label> <input type='text' title='Enter the Day, Month and year separated by spaces' class='form-control clickable uk-dob' style='' id='uk-dob-map' name='dob-map' value='00 00 0000' /><br/><label style='' for='uk-tob-map' class='control-label defaultedValueTob '>TOB (if known)</label> <input type='text' class='form-control clickable uk-tob' title='Enter the Hour and minute separated by a space' style='' id='uk-tob-map' name='tob' value='00 00' /></span> </div> <label for='uk-address' style='float:left;' class='control-label'>Location:</label><div class='col-sm-10' style='width:auto;'><input style='margin-left:4px;  width :363px' type='text' class='form-control clickable' id='uk-address'/></div></form> </div> <div class='modal-body requiresOnline' style='margin-top: 0px;'> <div class='form-horizontal' style='width: 550px'> <div class='form-group'> </div> <div id='uk-picker' style='width: 100%; height: 250px; margin-bottom:0px; '></div> </div> </div> <div class='modal-footer'> <button type='button' data-dismiss='modal' class='btn btn-primary requiresOnline' id='SaveChanges' >Save Changes + Get natal chart positions</button>  <button type='button' class='btn btn-default' data-dismiss='modal'>Close</button> </div> </div> </div>";

        var maphtmltar = " <div class='modal-dialog' style='background-color: rgb(202, 222, 230);'> <div class='modal-content' style='background-color: rgb(202, 222, 230);'> <div class='modal-header' style='background-color: rgb(202, 222, 230);'> <button type='button' class='close' data-dismiss='modal' aria-label='Close'><span aria-hidden='true'>&times;</span></button> <h4 class='modal-title'>Enter Your Birth Details</h4> </div> <div class='detailsForm' style='pointer-events: none; height: 130px; margin-top: 10px;' > <form id='f3' style='width:800px;' name='f3' class='n-t-small' action='index.pl' method='get'> <div class='entry-holder-tar' > <span style='white-space:nowrap;'> <label style='padding-right: 20px; width:180px;' for='uk-dob-map-tar' class=' control-label defaultedValueDob' >DOB</label> <input type='text' title='Enter the Day, Month and year separated by spaces' class='form-control clickable uk-dob-tar' style='' id='uk-dob-map-tar' name='dob-tar' value='00 00 0000' /><br/><label style='' for='uk-tob-map-tar' class='control-label defaultedValueTob '>TOB (if known)</label> <input type='text' class='form-control clickable uk-tob-tar' title='Enter the Hour and minute separated by a space' style='' id='uk-tob-map-tar' name='tob-tar' value='00 00' /></span> </div><label for='uk-address-tar' style='float:left;' class='control-label'>Location:</label> <div class='col-sm-10' style='width:auto;'><input style='margin-left: 4px;  width: 363px' type='text' class='form-control clickable' id='uk-address-tar'/></div> </form> </div> <div class='modal-body requiresOnline' style='margin-top: -20px;'> <div class='form-horizontal' style='width: 550px'> <div class='form-group'>  </div> <div id='uk-picker-tar' style='width: 100%; height: 250px; margin-bottom:0px; '></div> </div> </div> <div class='modal-footer'> <button type='button' data-dismiss='modal' class='btn btn-primary requiresOnline' id='SaveChangesTar' >Save Changes + Get natal chart positions</button>  <button type='button' class='btn btn-default' data-dismiss='modal'>Close</button> </div> </div> </div>";

        /*
         var maphtmltar="<div class='modal-dialog' style='background-color: rgb(202, 222, 230);'> <div class='modal-content' style='background-color: rgb(202, 222, 230);'> <div class='modal-header' style='background-color: rgb(202, 222, 230);'> <button type='button' class='close' data-dismiss='modal' aria-label='Close'><span aria-hidden='true'>&times;</span></button> <h4 class='modal-title'>Enter Your Birth Details</h4> </div> <div class='detailsForm' style='pointer-events: none; height: 130px; margin-top: 10px;' > <form id='f2' style='width:500px;' name='f2' class='n-t-small' action='index.pl' method='get'> <div class='entry-holder-tar' > <span style='white-space:nowrap;'> <label style='padding-right: 20px' for='uk-dob-map-tar' class='p-r-small col-sm-1 control-label defaultedValueDob' >DOB</label> <input type='text' title='Enter the Day, Month and year separated by spaces' class='form-control clickable uk-dob-tar' style='width: 150px; padding-right: 20px; margin-left: 50px; margin-top: -29px; float: left;' id='uk-dob-map-tar' name='dob-tar' value='00 00 0000' /> <label style='margin-left: 20px; padding-right: 50px;' for='uk-tob-map-tar' class='p-r-small col-sm-2 control-label defaultedValueTob '>TOB (if known)</label> <input type='text' class='form-control clickable uk-tob-tar' title='Enter the Hour and minute separated by a space' style='width: 110px; margin-left: 180px; margin-top: -30px; float: left;' id='uk-tob-map-tar' name='tob-tar' value='00 00' /></span> </div> </form> </div> <div class='modal-body requiresOnline' style='margin-top: -50px;'> <div class='form-horizontal' style='width: 550px'> <div class='form-group'> <label for='uk-address-tar' style='margin-left: -8px; margin-top: -50px;' class='col-sm-2 control-label'>Location:</label><div class='col-sm-10'><input style='margin-left: 80px; margin-top: -50px; width: 370px' type='text' class='form-control' id='uk-address-tar'/></div> </div> <div id='uk-picker-tar' style='width: 100%; height: 250px; margin-bottom:0px; '></div> </div> </div> <div class='modal-footer'> <button type='button' data-dismiss='modal' class='btn btn-primary requiresOnline' id='SaveChanges' >Save Changes + Get natal chart positions</button>  <button type='button' class='btn btn-default' data-dismiss='modal'>Close</button> </div> </div> </div>";
        */
        $('#uk-dialog').html(maphtml);
		
		

        $('#uk-dialog-tar').html(maphtmltar);

    }

    /* Assign the click  */
    $('#SaveChangesTar').on('click', function () {
        _f__SaveChangesGetChart("target", true);
    });

    /* Assign the click  */
    $('#SaveChanges').on('click', function () {
        _f__SaveChangesGetChart("subject", true);
    });

    var currentLat;
    var currentLon;
    var currentDob;
    var currentTob;

 
    currentLat = $('#uk-lat').prop('value');
    currentLon = $('#uk-lon').prop('value');
    currentDob = $('#uk-dob').prop('value');
    currentTob = $('#uk-tob').prop('value');



    $('#uk-picker').locationpicker({
        location: {
            latitude: currentLat,
            longitude: currentLon
        },
        radius: 0,
        inputBinding: {
            latitudeInput: $('#uk-lat'),
            longitudeInput: $('#uk-lon'),
            locationNameInput: $('#uk-address')
        },
        enableAutocomplete: true,
        onchanged: function (currentLocation) {

        }
    });
	
	
	
	$('#uk-picker').locationpicker('start');
	
	$('#uk-dialog,#uk-dialog-tar').on('hidden.bs.modal', function () {
   DOING_MAP_LOOKUP = false;
})
	
    $('#uk-dialog').on('shown.bs.modal', function () {
  DOING_MAP_LOOKUP = true;
        currentLat = $('#uk-lat').prop('value');
        currentLon = $('#uk-lon').prop('value')
        currentDob = $('#uk-dob').prop('value');
        currentTob = $('#uk-tob').prop('value');

        $('.uk-dob').prop('value', currentDob);
        $('.uk-tob').prop('value', currentTob);


        if (currentLat == 0 & currentLon == 0 & currentDob == '00 00 0000' & currentTob == '00 00') {
			
			 $('#uk-picker').locationpicker('location', {
            latitude: 51.4826,
            longitude: currentLon
        });
			
        }
		else{
			
			 $('#uk-picker').locationpicker('location', {
            latitude: currentLat,
            longitude: currentLon
        });
			
		}


       


        $('#uk-picker').locationpicker('autosize');


    });


    /*  Do the Target one  */

    currentLat = $('#uk-lat-tar').prop('value');
    currentLon = $('#uk-lon-tar').prop('value');

    $('#uk-picker-tar').locationpicker({
        location: {
            latitude: currentLat,
            longitude: currentLon
        },
        radius: 0,
        inputBinding: {
            latitudeInput: $('#uk-lat-tar'),
            longitudeInput: $('#uk-lon-tar'),
            locationNameInput: $('#uk-address-tar')
        },
        enableAutocomplete: true,
        onchanged: function (currentLocation) {

        }
    });
  


		$('#uk-picker-tar').locationpicker('start');
	
	
	$('#uk-dialog-tar').on('shown.bs.modal', function () {

	 DOING_MAP_LOOKUP = true;
	
        currentLat = $('#uk-lat-tar').prop('value');
        currentLon = $('#uk-lon-tar').prop('value')
        currentDob = $('#uk-dob-tar').prop('value');
        currentTob = $('#uk-tob-tar').prop('value');


        $('.uk-dob-tar').prop('value', currentDob);
        $('.uk-tob-tar').prop('value', currentTob);


        if (currentLat == 0 & currentLon == 0 & currentDob == '00 00 0000' & currentTob == '00 00') {
			  $('#uk-picker-tar').locationpicker('location', {
            latitude: 51.4826,  // Greenwich
            longitude: currentLon
        });
			
        }
		else
		{
			  $('#uk-picker-tar').locationpicker('location', {
            latitude: currentLat,
            longitude: currentLon
        });
			
		}

      
        $('#uk-picker-tar').locationpicker('autosize');
    });


    $('.hmsInputTypes, .bgTypes, .clearForm, .precessionFlag').selectmenu();
    $('.hmsInputTypes , .bgTypes, .clearForm, .precessionFlag').next().addClass("hmsInputTypesExtra");

    $(function () {
        $(".SignPrefix").selectmenu({
            change: function (event, ui) {
                var selected_value = ui.item.value;

                var TLA = SignNeatTitles[selected_value / 30];

                var theSignValueBox = $(this).next().next().next().next();

                theSignValueBox.prop('title', TLA);
                theSignValueBox.html(TLA.slice(0, 3));

                theSignValueBox.removeClass("unset-sign");

            }
        });
    });


    $("#selectionsHolder").dialog({
        autoOpen: false,
        width: 280,
        minWidth: 280,
        height: 645,
        minHeight: 605,
        modal: false,
        resizable: true,
        title: "Selections",
        dialogClass: "selectionsHolder",
        position: {
            my: "left bottom",
            at: "right top",
            of: "#logo"
        },
        show: {
            effect: "blind",
            duration: 300
        },
        hide: {
            effect: "blind",
            duration: 300
        },
        buttons: {
            /*  "Save Current Selection": function() {
			
			
        SaveCurrentSelection();
        }, */
            "Close": function () {
                $(this).dialog("close");
            }
        }

    });


    $('#selection_name').keyup(function () {

        if (typeof $(this).val() === "undefined")
            return;

        var valThis = $(this).val().toLowerCase();
        if (valThis == "") {
            $('#selections > li').show();
        } else {
            $('#selections > li').each(function () {
                var text = $(this).text().toLowerCase();
                (text.indexOf(valThis) >= 0) ? $(this).show(): $(this).hide();
            });
        };
    }).focus(function () {



    })

    ;


    /* 		//this allows us to bind events to the highly aloof JqueryUI select widget. Bleugh
    	$(function() {
                 $( "#selections" ).selectmenu({
    				 width: 360,
    				 height: 400,
    				 
                     change: function( event, ui ) {
    					 
    					 LoadChosenSelection(ui.item.value);

                     }
                 });
            });
    	 */


    $('.hmsbox') //leave submit buttons, records area, etc alone
        .bind('focus', function () {
            if ($(this).val() == "0 0 0") {
                $(this).val('');
            };
            $(this).addClass('greened');
        }) //more chaining = less searching
        .bind('blur', function () {
            _f__ConvertHMSToDecimal(this); /*  supposing you want to revert back on blur...*/
            $(this).removeClass('greened');
        })
        .bind('keyup', function () {

            $(this).val($(this).val().replace(/[^\d.]/g, ' ').replace(/ +/g, ' ').replace(/^ +/, '').replace(/^ /, ''));
            if ($(this).val() == '') {
                $(this).next().next().addClass("unset-sign");
            }
            _f__ConvertHMSToDecimal(this);

        });


    _f__clearCharts();

    $("#submit_btn").on("click", function () {
        _f__SaveChangesGetChart("subject", false)
    });
    $("#submit_btn-tar").on("click", function () {
        _f__SaveChangesGetChart("target", false)
    });
    _f__CheckCcCookieOptions(cc);


    $('.recordTitle').html('Activate xProfile mode');


    /* if( ! $('.xprofile-image').first().is(":visible")   )
    	{
    		$('.recordTitle').html('Activate chart mode');	
    	}
    	else
    	{
    		
    		
    	} */

}
function _f__populateBoxesFromDatastream(NameOfColumn) {

    if (!datastream.IsDefault()) {

        var passedData = datastream.split('¬');
        var boxes = $('.numVals-' + NameOfColumn);
        var simpleLocalIndex = 0;



 if (NameOfColumn == "subject") {
		SubjectsValues=[];
}else{
TargetsValues=[];
}


        jQuery.each(passedData, function (key, index) {
            if (key < 10 || key == 26 || key == 27) { // Sun -> Pluto, Asc., MC
                var splitted = index.split(',');

                boxes[simpleLocalIndex].value = (splitted[1]);
if (NameOfColumn == "subject") {
				SubjectsValues[SubjectsValues.length]=splitted[1];
}else
{TargetsValues[TargetsValues.length]=splitted[1];
}
                _f__ConvertDecimalToHMSSingle(boxes[simpleLocalIndex]);
                simpleLocalIndex++;
            }
        });


        var activeLonNode = "";
        var activeLatNode = "";
        var activeDobNode = "";
        var activeTobNode = "";
        var activeAccNode = "";
        if (NameOfColumn == "subject") {
           
            activeNameNode = $('#SubjectName');
            activeLonNode = $('#uk-lon');
            activeLatNode = $('#uk-lat');
            activeAccNode = $('#uk-acc');

        } else {
          
            activeNameNode = $('#TargetName');
            activeLonNode = $('#uk-lon-tar');
            activeLatNode = $('#uk-lat-tar');
            activeAccNode = $('#uk-acc-tar');
        }



/*
        if (!lon.IsDefault())
            activeLonNode.val(lat);
        if (!lat.IsDefault())
            activeLatNode.val(lon);
*/



        /*       
	   	$('#uk-dob').val(dob.replace('-',' ').replace('-',' '));
       	$('#uk-tob').val(tob.replace(':',' '));
*/
    }

    $('.error').hide();
	 activeNameNode.focus();

}

function _f__SetOptionsObjectFromWidgetState() {
    var hasError = 0;

    $.each(OptionsObject, function (index, value) {
        var radioButtons = $('#optionsHolder input:radio[name=' + index + ']');

        var selectedIndex = radioButtons.index(radioButtons.filter(':checked'));

        if (selectedIndex < 0) {
            hasError++;

            $("input:radio[name=" + index + "]:first").attr('checked', true);



            var radioButtons = $('#optionsHolder input:radio[name=' + index + ']');

            var selectedIndex = radioButtons.index(radioButtons.filter(':checked'));

            /* 	alert(selectedIndex);
            	
            	_f__warnDefaulted(index, "0"); */

            OptionsObject[index] = "0";
        } else {
            OptionsObject[index] = radioButtons.filter(':checked').prop('value');
        }


    });

    if (hasError > 0) {
        _f__SaveOptionsCookie(OptionsObject);
        _f__GetGlobalOptions(OptionsObject);
        /* _f__warnP("A recoverable error was encountered in the physical widget settings\n Please try again"); */
    }

    /*var statement = "";
	  $.each(OptionsObject, function  (index, value) {
	  statement += index + " is "+value+"\n";  
	  });
	alert (statement);*/
}

function _f__SetWidgetStateFromOptionsObject(reference) {

    $.each(reference,

        function (index, value) {
            $("input[name=" + index + "][value=" + value + "]").prop('checked', true);

            var indexCopy = index.toString();

            var forbiddenChars = new RegExp("[^a-zA-Z0-9]", 'g');
            if (forbiddenChars.test(index)) {
                indexCopy = index.replace(forbiddenChars, '');
            }


            if (!$.isNumeric(value)) {
                _f__warnDefaulted(index, value);
                value = "0"
            };
            if (parseFloat(value) < 0) {
                _f__warnDefaulted(index, value);
                value = "0"
            };

            _f__SetTovalue('OptionsObject.' + indexCopy, value);
            _f__SetTovalue(indexCopy, value);
        })
}

function _f__SetTovalue(varString, value) {
    eval(varString + " = " + value);
}

var allfunctions = [];

function _f__GetAllFunctions() {

    for (var i in window) {
        if ((typeof window[i]).toString() == "function") {
            allfunctions.push(window[i].name);
        }
    }

    allfunctions.sort(function (a, b) {
        return a.toUpperCase().localeCompare(b.toUpperCase());
    });

    var holder = [];

    allfunctions = $.grep(allfunctions, function (n, index) {
        return n.includes('_f__')
    });

    PrintArray(allfunctions);
}


function _f__DoFinalChartSteps (NameOfColumn){
	
	  	

        if (lat.IsDefault) {
            lat = activeLatNode.val();
        }
        if (lon.IsDefault) {
            lon = activeLonNode.val();
        };

        if (acc.IsDefault) {
            acc = 1 * (activeAccNode.is(':checked'))
        };

        if (lat.length == 0) {
            lat = "0";
            activeLatNode.focus();
            _f__warnP("Enter a decimal latitude value");
            return;
        };
        if (lon.length == 0) {
            lon = "0";
            activeLonNode.focus();
            _f__warnP("Enter a decimal longitude value");
            return;
        };

        /* 	var dayOB = splitted[0];
        var monOB = splitted[1];
        var yearOB = splitted[2];

        var minOB = splittedTime[1];
        var hourOB = splittedTime[0]; */

        var dataString = 'dob=' + dayOB + '/' + monOB + '/' + yearOB + '&tob=' + hourOB + ":" + minOB + '&lat=' + lat.toString().trim() + '&lon=' + lon.toString().trim() + '&acc=' + acc.toString();
        if (IsSubjectSave) {

            $('#uk-dob').val(dayOB + ' ' + monOB + ' ' + yearOB);
            $('#uk-tob').val(hourOB + ' ' + minOB);
        } else {
            $('#uk-dob-tar').val(dayOB + ' ' + monOB + ' ' + yearOB);
            $('#uk-tob-tar').val(hourOB + ' ' + minOB);
        }

        $.ajax({

            type: "POST",
            url: "./cgi-bin/chart.pl?",
            data: dataString,
            success: function (response) {
                $("#ScriptingInjectionPoint").html(response);

                /* alert( "Load was performed. dob:"+dob+" tob:"+tob+" lon:"+lon+" lat:"+lat+" response:"+response   ); */

                setTimeout(function () {
                    _f__populateBoxesFromDatastream(NameOfColumn);

					_f__LoadSubjectChartData(activeNameNode.val());

					GlobalCurrentSaveOptions[0] = activeLonNode.val();
					GlobalCurrentSaveOptions[1] = activeLatNode.val();
					GlobalCurrentSaveOptions[2] = activeDobNode.val() + " " + activeTobNode.val();
					GlobalCurrentSaveOptions[3] = activeAccNode.is(':checked') | 0;

					_f__doConditionalSave(NameOfColumn, activeNameNode.val(), GlobalCurrentSaveOptions);


                }, 150);
                return false;
            }
        });
    }

 var IsSubjectSave = false;

  var humDate = new Date();
    var timezone = "";
    var error = 0;

    var activeLonNode = "";
    var activeLatNode = "";
    var activeDobNode = "";
    var activeTobNode = "";
    var activeAccNode = "";
 
 
    var splitted = "";
    var splittedTime = "";

	   var timestamped =  ""; 
        var TimeZoneName = "";
        var TimeZoneId = "";
        var RawOffset = 0;
        var DstOffset = 0;
        var TotalOffset = 0;
	
	  var dayOB = "";
        var monOB = "";
        var yearOB = "";

        var minOB = "";
        var hourOB = "";
	
function _f__SaveChangesGetChart(NameOfColumn, useMap, passThisName) {

   if (typeof passThisName === 'undefined') { passThisName = 'default'; }
	
	IsSubjectSave = false;

     humDate = new Date();
     timezone = "";
     error = 0;

    activeLonNode = "";
    activeLatNode = "";
    activeDobNode = "";
    activeTobNode = "";
    activeAccNode = "";
 
 
   splitted = "";
   splittedTime = "";
	
	
    if (NameOfColumn.toLowerCase()  == "subject") {
        IsSubjectSave = true;

    }

	var activeNameNode

	if( IsSubjectSave){
			activeNameNode = $('#SubjectName')
        } else {
            activeNameNode = $('#TargetName')
        };
   
    if (IsSubjectSave) {
        if (useMap) {
            activeDobNode = $('#uk-dob-map')
        } else {
            activeDobNode = $('#uk-dob')
        };
        if (useMap) {
            activeTobNode = $('#uk-tob-map')
        } else {
            activeTobNode = $('#uk-tob')
        };
        activeLonNode = $('#uk-lon');
        activeLatNode = $('#uk-lat');
        activeAccNode = $('#uk-acc');
    } else {
        if (useMap) {
            activeDobNode = $('#uk-dob-map-tar')
        } else {
            activeDobNode = $('#uk-dob-tar')
        };
        if (useMap) {
            activeTobNode = $('#uk-tob-map-tar')
        } else {
            activeTobNode = $('#uk-tob-tar')
        };
        activeLonNode = $('#uk-lon-tar');
        activeLatNode = $('#uk-lat-tar');
        activeAccNode = $('#uk-acc-tar');
    }

    if (activeLatNode.val() == "") {
        error = 1;
        _f__warnClear();
        _f__warnP("Enter valid decimal latitude, e.g &nbsp; &nbsp; &nbsp; <strong>53.3</strong>");
        activeLatNode.focus();
        return;
    }
    if (activeLonNode.val() == "") {
        error = 1;
        _f__warnClear();
        _f__warnP("Enter valid decimal longitude, e.g &nbsp; &nbsp; &nbsp; <strong>-1.3</strong>");
        activeLonNode.focus();
        return;
    }

    

    splitted = activeDobNode.val().split(' ');
    splittedTime = activeTobNode.val().split(' ');
    splitted = $.grep(splitted, function (n) {
        return n.length > 0 || n
    });
    splittedTime = $.grep(splittedTime, function (n) {
        return n.length > 0 || n
    });
    if ((activeDobNode.val() == "00 00 0000") && (activeTobNode.val() == "00 00")) {
        error = 1;

	
        _f__warnP("ENTER VALID DOB/TOB (dd mm yyyy ) / (hh[ mm[ ss]]) for "+passThisName + " not "+activeDobNode.val() )

        return;
    }
    activeDobNode.val(activeDobNode.val());
    activeTobNode.val(activeTobNode.val());
    if (splittedTime.length != 2) {
        error = 1;
        _f__warnP("ENTER TIME OF BIRTH\n\nHH MM\n\nLike 13 00\n\nfor 1 PM for "+passThisName );
        return;
    }
    if (splitted.length != 3) {
        error = 1;
        _f__warnP("ENTER DATE OF BIRTH\n\nDD MM YYYY\n\nLike 11 09 1974\n\nfor 11th Sept. 1974 for "+passThisName );
        return;
    }
    if (error == 1) {
        _f__warnClear();
        _f__warnP("Please enter valid birth data for "+passThisName );
        return;
    }
    if (error == 0) {
         dayOB = splitted[0];
         monOB = splitted[1];
         yearOB = splitted[2];

         minOB = splittedTime[1];
         hourOB = splittedTime[0];

        OptionsObjectDateTime = new Date(Date.UTC(yearOB, (parseFloat(monOB) - 1), dayOB, hourOB, minOB));
        humDate = new Date(Date.UTC(yearOB, (parseFloat(monOB) - 1), dayOB, hourOB, minOB));

          timestamped = moment(humDate).unix(); 
         TimeZoneName = "";
         TimeZoneId = "";
         RawOffset = 0;
         DstOffset = 0;
         TotalOffset = 0;
	
		
		
		
             var url = "https://maps.googleapis.com/maps/api/timezone/json?location="+activeLatNode.val()+","+activeLonNode.val()+"&timestamp=" + timestamped ;
            $.ajax({
              url: url,
            }).done(function(response) {
 
		if (response.status == 'OK')
		{

        	TimeZoneName = response.timeZoneName.toString();
        	TimeZoneId = response.timeZoneId.toString();

        	if ($.isNumeric(parseFloat(response.rawOffset))){  OptionsObjectUtOffset=  parseFloat(response.rawOffset); RawOffset +=  parseFloat(response.rawOffset); }
        
        	if ($.isNumeric( parseFloat(response.dstOffset)) ) {OptionsObjectDstOffset=parseFloat(response.dstOffset); DstOffset +=  parseFloat(response.dstOffset);	} 

		}
		
		
        TotalOffset += RawOffset;
		TotalOffset -= DstOffset;
       humDate.setSeconds(humDate.getSeconds() + TotalOffset);
		
		
		
        	
 

		if (response.status == 'OK')
		{
			if (passThisName == 'default')
			  _f__infoP("BST offset: " +DstOffset/3600+ ",\n\n Total offset for that place/date is "+TotalOffset/3600+ " hour(s)\n\n", 2000); 
			else
			{_f__infoP ('importing '+passThisName, 2500 );	}	
		}
		else
		{
			if (passThisName == 'default')
			 _f__warnP("No timezone data for that location so offset assumed as 0 \n\n Total offset for that place/date is  0 secs for "+activeNameNode.val() , 1000); 
			else
			{	_f__infoP ('importing '+passThisName + ' with odd timezone', 2500 );}
		}

	
	
		setTimeout(function(){
 _f__DoFinalChartSteps(NameOfColumn);
 
 _f__populateRecordsList();
 
			}, (time_delay  * 0.5) );
		
        },0);
	}
}

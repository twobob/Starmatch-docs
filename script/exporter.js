var NEWARRAYLENGTH = 5;	// what does this control? old value 5

var exportDataList = [];
var exportedData = [];
var exportedLegibleDataString = "";
var exportString = "";
 /* <a href="../cgi-bin/exporter.pl"> RUN IT!</a> */

function _f__buildExporterData() {
    $('nav select option:selected').each(function (index, value) {
        exportDataList[exportDataList.length] = value.text; //.replace(' ', '-');
		console.log(exportDataList[exportDataList.length]);
        /* 	$(value).data(value.text.replace(' ','-'), value.value); */

		 exportedData[exportedData.length] = [value.text, value.value];

		
		
    })
	
	if (exportDataList.length == 0 ) 
	{
	
	   $('nav select option').each(function (index, value) {
        exportDataList[exportDataList.length] = value.text; //.replace(' ', '-');
        /* 	$(value).data(value.text.replace(' ','-'), value.value); */
	
	 exportedData[exportedData.length] = [value.text, value.value];
	
console.log(exportedData[exportedData.length-1]);
    })
	}
 	
	
	

    /* alert (exportedData); */

    /* 	quote Message-ID: <#C**************@TK2MSFTNGP12.phx.gbl>
    In JScript, variant string variables have the same limit as in VBScript,
    up to 2^31 characters.

    String *literals* (as in "this is a literal") have (IIRC) a limit ~2^10
    (1024) characters. */

    $.each(exportedData, function (index, value) {

	
        exportString += value.join('|') + delimeter;

    })

    exportString = exportString.slice(0, -1);

	
	
	
    /* alert(exportString); */
}
   
function _f__ajax_text() {
    var serverData = "";
    var input = exportString;
    $.ajax({
        url: "../cgi-bin/exporter.pl",
        type: 'POST',
        cache: false,
        data: {
            backup: input
        },
        success: function (serverData) {
            _f__handle_return(serverData);
        }
    });
}

var USEFUL_FILE_END_DELIMETER  = "\n";

function _f__create_flat(returnedString) {
	
	exportString = tempstring; //.replace(re, USEFUL_FILE_END_DELIMETER);;
	
}



var nameSelectionList = [];


function _f__handle_return(returnedString) {
	

	nameSelectionList = [];
	
	
	
	var updateCount = 0;
	
	
    var uncoded = returnedString.split(delimeter);

    var name = "";

	var nameList = "";
	
	
	
    $.each(uncoded, function (index, val) {

        var result = "";

		
		
		
	
	
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
		
		
		
		
		nameSelectionList[nameSelectionList.length] =  name;
		
		
		var addon = ", ";
		if(updateCount % 3 == 0)
		{
			addon = " ...<br />";
		}
		nameList += "["+updateCount+": "+ name+"] "+addon;
		
		_f__StorageSaveOperation(name,SubjectsValues,subjectProcessedData,myShakyOptionsArray   );

    })
	
	_f__updatedP("Updated "+updateCount+" Records for<br />"+nameList);
	
	
	
	
	setTimeout(function(){__f_SelectAndSave() }, 1);
	
	
	
}

function __f_SelectAndSave(){
	__f_SelectList(nameSelectionList); SaveCurrentSelectionWithName(fileInput.files[0].name.substr(0, fileInput.files[0].name.length-4));
	
	
}


function __f_SelectList(listToSelect){
	
	$(listToSelect).each( function(ind, val){     __f_SelectbyName(val)    }  )
	
}


function  __f_SelectbyName(nameToSelect)
{
	
	$('#person option').each( function(ind, val){  if (  $(val).text() == nameToSelect   ){    $(val).prop('selected', true)     }                        }  )
	
}




function _f__download(filename) {
	_f__buildExporterData();
	
  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(exportString));
  element.setAttribute('download', filename);

  element.style.display = 'none';
  document.body.appendChild(element);

  element.click();

  document.body.removeChild(element);
}


function _f__downloadLegible(filename) {
	_f__buildExporterData();
	
	
  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(exportString));
  element.setAttribute('download', filename);

  element.style.display = 'none';
  document.body.appendChild(element);

  element.click();

  document.body.removeChild(element);
}


function _f__createDownloadableFile(filename, passedData) {
  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(passedData));
  element.setAttribute('download', filename);

  element.style.display = 'none';
  document.body.appendChild(element);

  element.click();

  document.body.removeChild(element);
  
}


function _f__SendSaveRequest() {
    $.ajax({
        type: "POST",
        cache: false,
        url: "exporter.pl",
        processData: false,
        contentType: 'application/json',
        dataType: 'json',
        data: {
            backup: exportString
        }, // multiple data sent using ajax
        success: function (html) {

            /*   $('#add').val('data sent sent'); */

            alert(html);

            /*    $('#msg').html(html); */
        }
    });
    return false;

}

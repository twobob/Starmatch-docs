
var precision = 5;	// this is moved here so I can find it :( Will

function _f__setTextFactor(factorToSet)
{
	var refer = eval('NonRadioOptionsObject.'+factorToSet);
	var chosenInput = prompt("Please enter smam flag [0/1]",refer );
	
	if (chosenInput != null) {
		_f__SetTovalue('NonRadioOptionsObject.' + factorToSet, chosenInput);
		_f__SetTovalue(factorToSet, chosenInput);
	}
}

function _f__errorHighlight(e, type, icon) {
    if (!icon) {
        if (type === 'highlight') {
            icon = 'ui-icon-info';
        } else {
            icon = 'ui-icon-alert';
        }
    }
    return e.each(function  () {
        $(this).addClass('ui-widget');
        var alertHtml = '<div class="ui-state-' + type + ' ui-corner-all" style="padding:0.7em; background-color: #F2F2D0;">';
        alertHtml += '<span class="ui-icon ' + icon + '" style="  margin-top:4px;float:left;margin-right:.3em;"></span>';
        alertHtml += '<span id="out">';
        alertHtml += $(this).text() + "<br />";
        alertHtml += '</span>';
        alertHtml += '</div>';
        $(this).html(alertHtml);
    });
}

(function  ($) {
    $.fn.error = function  () {
        errorHighlight(this, 'error');
    };
})(jQuery);
(function  ($) {
    $.fn.highlight = function  () {
        errorHighlight(this, 'highlight');
    };
})(jQuery);


function guid() {
  function s4() {
    return Math.floor((1 + Math.random()) * 0x10000)
      .toString(16)
      .substring(1);
  }
  return s4() + s4() + '-' + s4() + '-' + s4() + '-' +
    s4() + '-' + s4() + s4() + s4();
}

function _f__flagP(arrayToPrint) {
	var realwidth =$(window).width();
	var modulus = 1;
	
	if( realwidth < 1080)
       modulus = 1;

    if (realwidth < 1500)
        modulus = 2;

	if( realwidth < 1890)
        modulus = 3;
	if ( realwidth > 1890)
        modulus = 4;
	
     var uuid = guid();

	 var textContent = "";
	var copyContent = "";
	
	 console.log(realwidth);
	 
	 var l = 0;

	$.each(arrayToPrint , function ( key, thing  ){
	copyContent += thing.toString()+" ";
	textContent += _f__padLeft(thing, 57, "·");
	
	if (l % modulus ==0 )  // We want an extra break for the title anyway
		{  	textContent += br;  
		copyContent += CRLF;
		}
	else
		{	textContent+= " "; 	
		copyContent += CRLF;
		}
	
	l++;
} );




var htmlthing =   textContent + '<br /><button class="copy-btn" data-clipboard-text="'+copyContent+'">copy</button>';


/* var htmlthing = "<button onclick='javascript:NonRadioOptionsObject.sdFactor=1 ;	_f__SetWidgetStateFromOptionsObject();_f__DoCompleteOptionsSaveThenLoadCycle();'>Set sdFactor</buttton>"; */



toastr.options = {
  "closeButton": false,
  "debug": false,
  "newestOnTop": true,
  "progressBar": false,
  "positionClass": "toast-bottom-center",
  "preventDuplicates": true,
  "onclick": null,
  "showDuration": "5000",
  "hideDuration": "1000",
  "timeOut": "0",
  "extendedTimeOut": "1000",
  "showEasing": "swing",
  "hideEasing": "linear",
  "showMethod": "fadeIn",
  "hideMethod": "fadeOut"
}

/* alert(OptionsObject.precessionFlag); */

toastr["error"](htmlthing, "Report!");

}



function _f__warnP(StringToPrint) {

var textContent = StringToPrint.toString();

toastr.options = {
  "closeButton": false,
  "debug": false,
  "newestOnTop": true,
  "progressBar": false,
  "positionClass": "toast-bottom-center",
  "preventDuplicates": true,
  "onclick": null,
  "showDuration": "5000",
  "hideDuration": "1000",
  "timeOut": "0",
  "extendedTimeOut": "1000",
  "showEasing": "swing",
  "hideEasing": "linear",
  "showMethod": "fadeIn",
  "hideMethod": "fadeOut"
}

toastr["error"](textContent, "Warning!");

}




function _f__updatedP(StringToPrint) {

var textContent = StringToPrint.toString();

toastr.options = {
  "closeButton": false,
  "debug": false,
  "newestOnTop": false,
  "progressBar": false,
  "positionClass": "toast-bottom-center",
  "preventDuplicates": false,
  "onclick": null,
  "showDuration": "5000",
  "hideDuration": "2000",
  "timeOut": "0",
  "extendedTimeOut": "1000",
  "showEasing": "swing",
  "hideEasing": "linear",
  "showMethod": "fadeIn",
  "hideMethod": "fadeOut"
}

toastr["info"](textContent, "UPDATED!");

}

function _f__debugBoxP(StringToPrint) {
  
  var textContent = StringToPrint.toString();
  
toastr.options = {
  "closeButton": false,
  "debug": false,
  "newestOnTop": true,
  "progressBar": false,
  "positionClass": "toast-bottom-center",
  "preventDuplicates": true,
  "onclick": null,
  "showDuration": "5000",
  "hideDuration": "1000",
  "timeOut": "0",
  "extendedTimeOut": "0",
  "showEasing": "swing",
  "hideEasing": "linear",
  "showMethod": "fadeIn",
  "hideMethod": "fadeOut"
}

toastr["info"](textContent,"Debug");

}

function _f__infoP(StringToPrint, timeToClear) {
function _f__isNumeric(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}

var textContent = "";
if(Object.prototype.toString.call( StringToPrint ).includes('Array')) {
$.each(StringToPrint,function  (index, value) {
		textContent += value+", ";
} )
}
else
{
textContent = StringToPrint.toString();
}
if (typeof timeToClear === "undefined" )
{
timeToClear=0
}

toastr.options = {
  "closeButton": false,
  "debug": false,
  "newestOnTop": true,
  "progressBar": false,
  "positionClass": "toast-bottom-center",
  "preventDuplicates": true,
  "onclick": null,
  "showDuration": "0",
  "hideDuration": "1000",
  "extendedTimeOut": timeToClear,
  "showEasing": "swing",
  "hideEasing": "linear",
  "showMethod": "fadeIn",
  "hideMethod": "fadeOut"
}

toastr.options.extendedTimeOut = timeToClear; //1000;
    toastr.options.timeOut = timeToClear;

toastr["success"](textContent, "Info:");


}


 function _f__getLastToast(){
            return $toastlast;
        }
	
var repaintRequired = false;
var myTimeout;

function handleDebugRepaint(){
	
	repaintRequired = false
	$('#debugArea').show();
}	

function debugP(StringToPrint) {
	
	$('#debugArea').hide();
	
 	$('#debugArea').html($('#debugArea').html() +   StringToPrint.toString() + CRLF + "<br />");
	
	if (!repaintRequired)
	{
	 myTimeout = setTimeout(handleDebugRepaint, 200);
	}
	
}
		

function _f__debugClear() {
    $("#debugArea").hide();
    document.getElementById('debugArea').innerHTML = "";
	 $("#debugArea").show();
	toastr.clear();
}

function _f__warnClear() {

    $("#warningArea").height("0px");
    document.getElementById('warningArea').innerHTML = "";
	toastr.clear();
}


if (!Date.now) {
    Date.now = function  () {
        return new Date().getTime();
    }
}

function _f__autofill() {
    arrPositions = PlanetNeatTitles;
    for (i = 0; i < 12; i++) {
        document.f1.label[i].value = arrPositions[i];
        document.f1.numvalue1[i].value = parseFloat(Math.random() * 360).toFixed(precision);
        document.f1.numvalue2[i].value = parseFloat(Math.random() * 360).toFixed(precision);
    }
}

function _f__copy(o) {
    var out, v, key;
    out = Array.isArray(o) ? [] : {};
    for (key in o) {
        v = o[key];
        out[key] = (typeof v === "object") ? copy(v) : v;
    }
    return out;
}

function _f__encode(c) {
    var x = 'charCodeAt',
        b, e = {},
        f = c.split(""),
        d = [],
        a = f[0],
        g = 256;
    for (b = 1; b < f.length; b++) c = f[b], null != e[a + c] ? a += c : (d.push(1 < a.length ? e[a] : a[x](0)), e[a + c] = g, g++, a = c);
    d.push(1 < a.length ? e[a] : a[x](0));
    for (b = 0; b < d.length; b++) d[b] = String.fromCharCode(d[b]);
    return d.join("")
};

function _f__decode(b) {
    var a, e = {};
    var d = b.split("");
	var f = d[0];
        var c = d[0];
        var g = [c];
        var h = 256;
		var o = 256;
    for (b = 1; b < d.length; b++) a = d[b].charCodeAt(0), a = h > a ? d[b] : e[a] ? e[a] : f + c, g.push(a), c = a.charAt(0), e[o] = f + c, o++, f = a;
    return g.join("")
};

function _f__CleanupNaN(valueRecord) {
    if (valueRecord.val() == 'NaN') {
        valueRecord.val('0');
    }
} 
    


function _f__TogglexProfile() {

if( ! $('.xprofile-image').first().is(":visible")   )
	{
		
		$('.recordTitle').html('Activate chart mode');
		
		$('.hideInOneToOne').attr('style','display:block !important');
		$('.hideInxProfile').attr('style','display:none !important');
			
	
		
	}
	else
	{
			$('.recordTitle').html('Activate xProfile mode');
		
			$('.hideInxProfile').attr('style','display:block !important');
				$('.hideInOneToOne').attr('style','display:none !important');
				
				
				
		
	}


	$(this).blur();
    $('.hideableRecords').toggleClass("hiddenThing");


    $('.main-logo').toggleClass("smallerLogo");
    $('.mainForm').toggleClass("shrinkForm");
    $('#resizableOptions').toggleClass("resizableOptions");
	
	
	
	
	$('.TargetStuffRow').show();

}

	
function _f__toggleLookups(){
	
	debugP("_f__toggleLookups is deprecated");
	
	/* if(!  $('.xprofile-image').first().is(":visible")   )
	{
	$('.TargetStuffRow').toggle();

	}
	
	$('.details_form, .hiddenWhenMeta').toggle();
	

	
toggleLinkText(); */

}


function toggleLinkText(){
	
   	$('.targetArea').prop('disabled', !$('.targetArea').prop('disabled')    );
	
		
	if(  $('.xprofile-image').first().is(":visible")   )
	{
		$('.toggleMessageLink').html('');
		
	}
	else
	{
		
		$('.toggleMessageLink').html('Close this area to perform analyisis.');
		
	}
}	
	
function _f__ConvertHMSToDecimal(idOfThing) {
    var getthevalue = $(idOfThing).attr('id');

    var thehmsValueBox = $("#" + getthevalue);

    var i = thehmsValueBox.val().trim().split(" "),
        h = i[0],
        n = i[1],
        s = i[2];

    if (typeof (h) === "undefined") {

        h = n = s = 0;
 
        thehmsValueBox.val('0 0 0');
		thehmsValueBox.blur();
		thehmsValueBox.focus();
		thehmsValueBox.next().next().addClass("unset-sign");
    /*     showInfo(); */
	/* 	warnP(' &deg; must be in the range of 0 to <360. [undefined] '); */
        return;
    }

    if (h.length == 0) {
        thehmsValueBox.val('0 0 0');
		thehmsValueBox.blur();
		thehmsValueBox.focus();
		thehmsValueBox.next().next().addClass("unset-sign");
       /*  showInfo(); */

	/* 	warnP(' &deg; must be in the range of 0 to <360.  '); */
        return;

    }



    if (!$.isNumeric(h)) {

        warnClear();
       /*  warnP(' &deg; must be in the range of 0 to <360. [numeric] '); */
        thehmsValueBox.val('0 0 0');
		thehmsValueBox.next().next().addClass("unset-sign");
		thehmsValueBox.blur();
		thehmsValueBox.focus();
        return;
    }

	
	  absdlat = Math.abs(h);
	
    if (absdlat >= 360) {
        warnClear();
      /*     warnP(' &deg; must be in the range of 0 to <360. '); */
        thehmsValueBox.val('0 0 0');
		thehmsValueBox.next().next().addClass("unset-sign");
		thehmsValueBox.blur();
		thehmsValueBox.focus();
        return;
    }


    var selectedSignValue = 0;

    if (parseFloat(h) >= 30) //  WE WANT TO SET THE SIGN!

    {

        selectedSignValue = parseFloat(h);

    } else if (thehmsValueBox.prev().prev().val() > 0) {

	
	
        selectedSignValue = parseFloat(thehmsValueBox.prev().prev().val()) + (parseFloat(h) % 30);

    }
	else 
	{
		 selectedSignValue = parseFloat(h);
		
	}
	
	
	
    var temp = selectedSignValue;

    thehmsValueBox.prop('title', SignNeatTitles[Math.floor(selectedSignValue / 30)]);

	
	var TLA = SignNeatTitles[Math.floor(selectedSignValue / 30)];
		

	thehmsValueBox.next().next().prop('title', TLA);
	thehmsValueBox.next().next().html(TLA.slice(0,3));

	thehmsValueBox.next().next().removeClass("unset-sign");

	
    if (absdlat >= 30) {
        absdlat = absdlat % 30;
    }

    absdlat = Math.abs(absdlat * 1000000.); //integer
	
    if (typeof (n) === "undefined" || !$.isNumeric(n)) {
	
	
	/* warnClear();
	warnP("Enter &deg;Minutes"); */
		
	$(this).focus();
	return;
    /*     n = "0";
        s = "0"; */

    } else if (typeof (s) === "undefined" || !$.isNumeric(s)) {

      /*   s = "0"; */
    
		/* 
		warnClear();
	warnP("Enter &deg;Seconds"); */
$(this).focus();
	return;
    }
    


    s = Math.floor(parseFloat(s) + 0.5);
    thehmsValueBox.val(h + ' ' + n + ' ' + s);
	
	
    absmlat = Math.abs(n * 1000000.); //integer

    if (absmlat >= (60 * 1000000.) || !$.isNumeric(n)) {
    /*     warnClear();
        warnP(' &deg;mins must be in the range of 0 to <60. '); */
        thehmsValueBox.val(h + ' 0 0');
		
    }


    absslat = Math.abs(s * 1000000.); // Note: kept as big integer for now, even if submitted as decimal

    if (absslat >= (60 * 1000000.) || !$.isNumeric(s)) {
     /*    warnClear();
        warnP(' &deg;seconds must be in the range of 0 to <60. '); */
        thehmsValueBox.val(h + ' ' + n + ' 0');
        return;
    }

    absdlat = Math.abs(temp * 1000000.); //integer

    var FinalConvertedValue = Math.round(absdlat + (absmlat / 60.) + (absslat / 3600.)) * latsign / 1000000;


    if (absdlat >= 30 * 1000000.) {


        var signIndexer = Math.floor(selectedSignValue / 30);
        
		
		var TLA = SignNeatTitles[signIndexer];
		
		
		
		
		thehmsValueBox.next().next().prop('title', TLA);
	thehmsValueBox.next().next().html(TLA.slice(0,3));

	thehmsValueBox.next().next().removeClass("unset-sign");
	
		
		
		thehmsValueBox.prev().prev().prop('selectedIndex', signIndexer);


        thehmsValueBox.val(h % 30 + ' ' + n + ' ' + s);

    }


    thehmsValueBox.next().val(FinalConvertedValue);

}


function _f__ConvertDecimalToHMSSingle(idOfThing) {
    var stretch = 10000;

    var getthevalue = $(idOfThing).attr('id');
    var theDecimalValueBox = $("#" + getthevalue);


    if (!$.isNumeric(theDecimalValueBox.val())) {
        theDecimalValueBox.val(0);

        return;
    }

   var latAbs = Math.abs(Math.round(theDecimalValueBox.val() * stretch));
    var hours = (Math.floor(latAbs / stretch) * 1);


    theDecimalValueBox.prev().prev().prev().prop('selectedIndex', Math.floor(hours / 30));

    var mins = Math.floor(((latAbs / stretch) - Math.floor(latAbs / stretch)) * 60);
    var secs = (Math.floor(((((latAbs / stretch) - Math.floor(latAbs / stretch)) * 60) - Math.floor(((latAbs / stretch) - Math.floor(latAbs / stretch)) * 60)) * stretch) * 60 / stretch);
 

	if (hours + mins + secs > 0)
	{
	var TLA = SignNeatTitles[Math.floor(hours / 30)];	
		
	theDecimalValueBox.next().prop('title', TLA);
	theDecimalValueBox.next().html(TLA.slice(0,3));
		
	theDecimalValueBox.next().removeClass("unset-sign");
	}
	
	theDecimalValueBox.prev().val(hours % 30 + " " + mins + " " + Math.floor(secs + 0.5));

    theDecimalValueBox.prev().prop('title', SignNeatTitles[Math.floor(hours / 30)]);
}


function _f__showInfo() {

    warnClear();

}
	
function _f__AutoSetupTheTabbingOrder() {
    var tabindex = 1;
    $('.hmsbox,option').each(function  () {
        if (this.type != "hidden") {
            var $input = $(this);
            $input.attr("tabindex", tabindex);
            tabindex++;
        }
    });
    $('.noPeriod').on('input', function  () {
        var chars = $(this).val();
        if (chars.indexOf('.') !== -1) {
            chars = chars.replace('.', '');
            $(this).val(chars);
        }
    });


}


var emptySearchFilter = function(){ 
_f__handleKeyPress($('#filterTxt').val(), '');
$('#filterTxt').val('');
}


function _f__InitSearching() {
    var prevTxt = null;
    var txt = document.getElementsByName('filterTxt');
    if (txt != null) {
        txt[0].onkeyup = function  (event) {
            var e = event || window.event;
            var curTxt = txt[0].value;
            _f__handleKeyPress(prevTxt, curTxt);
            prevTxt = curTxt;
            return true;
        }
    }
}


function _f__handleKeyPress(oldVal, newVal) {
    var components = document.getElementsByName('person');
    var select = components[0];
    if (originalentries === null) {
        originalentries = new Array();
        for (c = 0; c < select.children.length; c++) {
            originalentries.push(select.children[c]);
        }
    }
    if (oldVal !== null && (newVal.length < oldVal.length)) {
        for (c = 0; c < originalentries.length; c++) {
            select.add(originalentries[c]);
        }
    }
    var parts = newVal.split(' ');
    var toremove = new Array();
    for (i = 0; i < select.children.length; i++) {
        var entry = select.children[i];
        var match = true;
        var entryTxt = entry.text;
        for (p = 0; p < parts.length; p++) {
            var part = parts[p].toUpperCase();
            if (entryTxt.toUpperCase().lastIndexOf(part) < 0) {
                match = false;
                break;
            }
        }
        if (match == false) {
            toremove.push(entry);
        }
    }
    if (toremove != null) {
        for (t = 0; t < toremove.length; t++) {
            var entryTxt = toremove[t].text;
            select.removeChild(toremove[t]);
        }
    }
	
	var TotalContentCount = $('#person option').length;
	
	debugP(TotalContentCount +" of " +  originalentries.length + " records");
	
	
}
 

 
function _f__ShowAdminAreas() {

	
	$( "#optionsHolder" ).dialog( "open" );
}

function _f_ShowExporter(){
	$( "#exporterHolder" ).dialog('open');
	
}




function _f__BindKeys() {

    Mousetrap.bind(['alt+shift+o'], function  (e) {
        _f__ShowAdminAreas();
        return false;
    });

    Mousetrap.bind(['shift+esc'], function  (e) {
		
        _f__debugClear();
_f__infoP("cleared", 1000); //		toastr.clear();
        return false;
    });
	
  Mousetrap.bind(['alt+shift+x'], function  (e) {
	
	$( "#exporterHolder" ).dialog('open');
	
 });
}

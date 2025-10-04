google.load("visualization", "1.1", {packages:["corechart", "table"]}); 
function _f__SetAnalysisMethod(num){
	method = parseInt(num);
	_f__ClickedxProfile();
}



function _f__drawChart(title, arrayToDraw, idOfDivToUse) {
	var data = "";
	var options = "";
	var seriesBarColor = "#518EAB";

	if(idOfDivToUse.includes("2"))
		seriesBarColor = "#71BEDB";

	if(!$("#chartsHolder").dialog( "isOpen" ))
		$("#chartsHolder").dialog( "open" );

	if ( OptionsObject.chartType==1 )
	{
		var stats = { mean: 0, sd: 0, sError: 0, skew: 0 };
		compositeThemeValues ( arrayToDraw );
		for ( n = 0 ; n < arrayToDraw.length; n++ )
			arrayToDraw[n] = avgThemeVal[n];

		statAnalyse ( arrayToDraw, 0, arrayToDraw.length, stats, 0 );	// pop. std. dev.
/*
debugP("output");
debugP("s mean "+pNum(stats.mean, 3));
debugP("s sd "+pNum(stats.sd, 3));
*/
		var lo = stats.mean-(stats.sd*sigmaFactor);
		var hi = stats.mean+(stats.sd*sigmaFactor);
		var mean = stats.mean ;
/*
debugP("output.js "+title+" type 1 statistics:");
debugP("output.js type 1 array:");
for ( n= 0 ; n < 12; n++ )
debugP(pNum(arrayToDraw[n], precision));
debugP("mean: "+pNum(mean, precision));
debugP("sigma: "+pNum(stats.sd, precision));
debugP("limits (+/ 1 sigma): "+pNum(lo, precision)+", "+pNum(hi, precision));
debugP("std error on mean: s = "+pNum(sStats.sError, precision)+" t = "+pNum(tStats.sError, precision));
debugP("");
*/

		data = google.visualization.arrayToDataTable(
        [
            ['theme', 'amount', {id:'i0', type: 'number', label: 'Mean', role: 'interval'}, {id:'i1',type: 'number', label: 'hi', role: 'interval', color: 'yellow'}, {id:'i2',type: 'number', label: 'lo', role: 'interval'},{id:'i0',type: 'number', label: '', role: 'interval'}],
            ['1', arrayToDraw[0],mean, lo, hi ,mean  ],
            ['2', arrayToDraw[1],mean, lo, hi ,mean  ],
            ['3', arrayToDraw[2],mean, lo, hi ,mean  ],
            ['4', arrayToDraw[3],mean, lo, hi   ,mean],
            ['5', arrayToDraw[4],mean, lo, hi   ,mean],
            ['6', arrayToDraw[5],mean, lo, hi   ,mean],
            ['7', arrayToDraw[6],mean, lo, hi   ,mean],
            ['8', arrayToDraw[7],mean, lo, hi   ,mean],
            ['9', arrayToDraw[8],mean, lo, hi   ,mean],
            ['10', arrayToDraw[9],mean, lo, hi  ,mean],
            ['11', arrayToDraw[10],mean, lo, hi ,mean],
            ['12', arrayToDraw[11],mean, lo, hi ,mean]
        ]);

		options = {
		width: 360,
		height: 200,
	
	interval: {
            'i0': { 'style':'line', 'color':'#D3362D', 'lineWidth': 0.5 },
            'i1': { 'style':'line', 'color':'#F1CA3A', 'lineWidth': 1 },
            'i2': { 'style':'line', 'color':'#5F9654', 'lineWidth': 2 },
        },

	
        orientation: 'horizontal',
		bar: {
		groupWidth: "99%"
        },
		seriesType: 'bars',
		intervals: {
			shortBarWidth: '5',
			barWidth: '5',
			pointSize : 0
		},
		
        series: {
			0: {type: 'bars',
				color: seriesBarColor
			},
			3: {
				type: 'line'
			},
			4: {
				type: 'line'
			}
		},
        enableInteractivity: true,
        backgroundColor: '#D3E0E8',
         title: title,
        legend: 'none',
        chartArea: {
			width: '100%',
			height: '80%'
		},
        hAxis: {
			showTextEvery: 1,
			title: '',
			minValue: 0
		},
        vAxis: {
			showTextEvery: 1,
			gridlines: {
				color: 'transparent'
			},
			title: ''
		}
	};
}
else  // Chart 0 and any chart not chart 1
{
/*
debugP("output: type 1 array");
for ( n= 0 ; n < 12; n++ )
debugP(pNum(arrayToDraw[n], precision));
*/
	normaliseValues ( arrayToDraw, arrayToDraw.length );

	var vCount =6;

     data = google.visualization.arrayToDataTable(
        [
            ['theme', 'amount'],
            ['1', arrayToDraw[0]],
            ['2', arrayToDraw[1]],
            ['3', arrayToDraw[2]],
            ['4', arrayToDraw[3]],
            ['5', arrayToDraw[4]],
            ['6', arrayToDraw[5]],
            ['7', arrayToDraw[6]],
            ['8', arrayToDraw[7]],
            ['9', arrayToDraw[8]],
            ['10', arrayToDraw[9]],
            ['11', arrayToDraw[10]],
            ['12', arrayToDraw[11]]
        ]);

		
    options = {
		 width: 360,
		 height: 200,
	
        orientation: 'horizontal',
        bar: {
            groupWidth: "99%"
        },
        series: {
            0: {
                color: seriesBarColor
            }
        },
        enableInteractivity: false,
        backgroundColor: '#D3E0E8',
        title: title,
        legend: 'none',
        chartArea: {
           width: '100%',
            height: '80%'
        },
        hAxis: {
            showTextEvery: 1,
            title: '',
            minValue: 0
        },
        vAxis: {
            showTextEvery: 1,
            gridlines: {
color: "red",
count: vCount
            },
            title: ''
        }
    };

		

}




	
    var chart = new google.visualization.BarChart(document.getElementById(idOfDivToUse));

	
    chart.draw(data, options);


	
}

function _f__makeBarColour(barColour, TOBaccurate){ 
var blueTint = 0.655;
var redTint = 0.314;
var greenTint = 0.557;

if (method>1 && method<4)
{
	if (barColour<0 && barColour>-4)
	{
		redTint = 1-(Math.abs(barColour+1)/4);
		greenTint = 0;
		blueTint = 0;
	}
	if (barColour>=0 && barColour<3)
	{
		greenTint = 1-(barColour/4);
		blueTint = 0;
		redTint = 0;
	}
}

if ( (method == 1) && (TOBaccurate == 'false') )	// recolour entries with flagged inaccurate bt
{	// a few shades lighter than the default
	redTint = 0.44;
	greenTint = 0.78;
	blueTint = 0.917;
}

redTint = Math.floor(redTint*255);
greenTint = Math.floor(greenTint*255);
blueTint = Math.floor(blueTint*255);
return 'rgb('+redTint+','+greenTint+','+blueTint+')';
}


function _f__drawxProfileChart(title, arrayToDraw, idOfDivToUse, stats) {

debugP("output.js _f__drawxProfileChart()");
debugP("Fit statistics:");
var mean;
mean = stats.mean;
var lo = stats.mean-(stats.sd*0.5);
var hi = stats.mean+(stats.sd*0.5);
debugP("mean "+pNum(stats.mean, precision)+" std. dev. "+pNum(stats.sd, precision)+" limits: lo "+pNum(lo, precision)+" hi "+pNum(hi, precision));

var mean = 0;
lo = -stats.sd*0.5;
hi = stats.sd*0.5;
debugP("adjusted mean "+pNum(mean, precision)+" std. dev. "+pNum((stats.sd*sigmaFactor), precision)+" limits: lo "+pNum(lo, precision)+" hi "+pNum(hi, precision));

	if(!$("#xchartsHolder").dialog( "isOpen" ))
		$("#xchartsHolder").dialog( "open" );

	$('.xchartsHolder').show();
	$('.chartsHolder').hide();
	var data = new google.visualization.DataTable();
    data.addColumn('string', 'Name');
    data.addColumn('number', 'Weight');
    data.addColumn({type: 'string', role: 'style'});
	data.addColumn({id:'i0', type: 'number', label: 'mean', role: 'interval'});
	data.addColumn({id:'i1', type: 'number', label: 'lo', role: 'interval'});
	data.addColumn({id:'i2', type: 'number', label: 'hi', role: 'interval'});
	data.addColumn({id:'i0', type: 'number', label: '', role: 'interval'}); // fudge...
	data.addColumn({type: 'string', role: 'annotation'});
	data.addColumn({type: 'string', role: 'annotationText'});

	$.each(arrayToDraw, function  (key, value) {
var expressed = _f__makeBarColour(value[2], value[3]);
	

data.addRow([value[0],value[1], expressed  ,    mean, lo, hi ,mean,       value[0].toString(), value[1].toFixed(3) +': '+value[0].toString()] );
     });

    var options = {
		tooltip: {
		isHtml:true
   
  },
		annotations: { alwaysOutside: false},
		
		interval: {
            'i0': { 'style':'line', 'color':'#D3362D', 'lineWidth': 0.5 },
            'i1': { 'style':'line', 'color':'#F1CA3A', 'lineWidth': 1 },
            'i2': { 'style':'line', 'color':'#5F9654', 'lineWidth': 2 },
        },
		intervals: {
			shortBarWidth: '5',
			barWidth: '5',
			pointSize : 0
		},
		'allowHtml': true,
		 width: '90%',
		 height:  '100%',
        orientation: 'vertical',
        bar: {
            groupWidth: "99%"
        },
        enableInteractivity: true,
        backgroundColor: '#D3E0E8',
        title: title,
        legend: 'none',
        chartArea: {
            width: '80%',
            height: '80%'
        },
      
        vAxis: {
       textPosition: 'none',
            title: ''
		
        },
        hAxis: {
			format: 'short'
          
		}
    };
    var chart = new google.visualization.BarChart(document.getElementById(idOfDivToUse));
    chart.draw(data, options);
/* 	 */$('.xchartsHolder').show();
}

function _f__drawWeightingChart(title, nameArray, arrayToDraw, idOfDivToUse) {

	
	
	  var data = google.visualization.arrayToDataTable(
        [
            ['Weighting', nameArray[0], nameArray[1], ],
            ['', arrayToDraw[0], arrayToDraw[1]]
        ]);

	var bg = '#d3e0e8';
	
	 var options = {
         width: 350,
	      height: '50',
        legend:  {position: 'none'},
         bar: { groupWidth: '100%' },
		 enableInteractivity: true,
        isStacked: true,
		  backgroundColor: '#d3e0e8',
		   colors:['#518EAB','#71bEdB'],
		  hAxis: {
			  	 baselineColor: bg,
         gridlineColor: bg,
         textPosition: 'none',
            showTextEvery: 1,
            title: ''
        },
        vAxis: {
			 baselineColor: bg,
         gridlineColor: bg,
         textPosition: 'none',
            showTextEvery: 1,
            gridlines: {
                color: 'transparent'
            },
            title: ''
        } 
    
		  
      };
	
   /*  var options = {
        orientation: 'vertical',
        bar: {
            groupWidth: "99%"
        },
        series: {
            0: {
                color: '#518EAB'
            }
        },
        enableInteractivity: true,
        backgroundColor: '#D3E0E8',
		isStacked: true,
		isStacked: 'percent',
        title: title,
        legend: 'none',
        chartArea: {
            width: '100%',
            height: '80%'
        },
        hAxis: {
            showTextEvery: 1,
            title: ''
        },
        vAxis: {
            showTextEvery: 1,
            gridlines: {
                color: 'transparent'
            },
            title: ''
        }
    }; */

    var chart = new google.visualization.BarChart(document.getElementById(idOfDivToUse));

    chart.draw(data, options);
}

function _f__GenerateChart() {
      if(!_f__RecordEntriesAreNoneZero(SubjectsValues))
	  {
		  
		  _f__infoP("Load a subject with sensible positions, Thank you");
		  return;
	  }
	    if(!_f__RecordEntriesAreNoneZero(TargetsValues))
	  {
		  
		  _f__infoP("Load a target with sensible positions, Thank you");
		  return;
	  }


if(!$("#chartsHolder").dialog( "isOpen" ))
	{
		
		$("#chartsHolder").dialog( "open" );
		
		
	}


	$('.chartsHolder').show();

	$('.xchartsHolder').hide();


	
    _f__warnClear();
	

	_f__DoCompleteOptionsSaveThenLoadCycle();

	
	_f__LoadBothChartsData( $('#SubjectName').val(), $('#TargetName').val()  );
   

    _f__drawChart($('#SubjectName').val(), SubjectProcessedData, 'chart_div_1');
	
    _f__drawChart($('#TargetName').val(), TargetsProcessedData, 'chart_div_2');

	

	
var subReporter = "Chart for " +$('#SubjectName').val();
	
var	reporter = [];
var sLength = SubjectProcessedData.length;
reporter[reporter.length] =subReporter;	


subReporter = "";
    for (l = 0; l < sLength; l++) {
		var sliced = SubjectProcessedData[l].toString().slice(0,6);

	if ( sliced.length>1 )
	{
		subReporter +=  _f__pad(sliced, 6, "0");
		if ( l != sLength-1 )
			subReporter += ",";
	}
	else
	{
		subReporter += sliced;
		if ( l != sLength-1 )
			subReporter += ",";
	}
}
reporter[reporter.length] =subReporter;

var tLength = TargetsProcessedData.length;
subReporter = "Chart for " +$('#TargetName').val();
reporter[reporter.length] =subReporter;	
subReporter = CRLF;
	
    for (l = 0; l < tLength; l++) {
		
	    
	    
		var sliced = TargetsProcessedData[l].toString().slice(0,6);
		
	        if ( sliced.length>1 )
		{
			subReporter +=  _f__pad(sliced, 6, "0");
			if ( l != tLength-1 )
				subReporter += ",";
		}

	        else
        	{
        	        subReporter += sliced;
			if ( l != tLength-1 )
	                	subReporter += ",";
        	}
    }

	reporter[reporter.length] =subReporter;
	

	lastFullChartReport = reporter;

	
/*
debugP("subject");
debugP(SubjectProcessedData);
debugP("target");
debugP(TargetsProcessedData);
*/

	relStrength(SubjectProcessedData,TargetsProcessedData );
var charWeighting = relStrength(SubjectProcessedData,TargetsProcessedData );

	_f__drawWeightingChart("Weighting", [$('#SubjectName').val(),$('#TargetName').val() ], charWeighting , 'chart_div_3' );
}

 function add(a, b) {   return a + b;  }
 function relStrength (SubjectData,TargetsData )
	{
		var Subtot =  SubjectData.reduce(add, 0), TargTot =  TargetsData.reduce(add, 0), TotTot = Subtot + TargTot;
		return [ (1.0 / TotTot) * Subtot , (1.0 / TotTot) * TargTot  ];
	}

function _f__clearCharts() {

	$('.chartsHolder').hide();
	$('.xchartsHolder').hide();
	
}

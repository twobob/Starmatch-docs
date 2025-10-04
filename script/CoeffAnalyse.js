var y0, y1;


function _f__coeffAnalysis() {	// called by global.js _f__InitPage()
	var v, x, y, y0, y1, y2;
	var curveData = { aQ: 0, bQ: 0,  cQ: 0, m0: 0, c0: 0, m1: 0, c1: 0, centre: 0, distL:0, distR: 0 };
	var ixL, ixR, xOffset;
	var i, j, n, m, numRecords;
	var numPrimaryPeaks = [0,0,0,0,0,0,0,0,0,0,0,0];
	var numPrimaryTroughs = [0,0,0,0,0,0,0,0,0,0,0,0];

	i = 0;
	j = 0;
	n = 0;
	m = 0;
	numRecords = curveDataSet[0];
debugP("CurveAnalyse: records: "+numRecords);	// number of subjects' graphs
	var dataSize, yval;
	while ( i < numRecords )
	{
		var numFunctions = curveDataSet[j+2];
dataSize = curveDataSet[j+3];
		var peak = 0;
		var peakLoc = -1;	// invalid
		var trough = 0;
		var troughLoc = -1;	// invalid
		for ( n = 0; n < numFunctions; n++ )
		{
			m = n+j+4;
			curveData = curveDataSet[m];
			yval = ordinate ( curveData.centre, curveData, dataSize )	// x is either peak or trough location
if ( curveData.aQ < 0 )	
{
			if ( yval > peak )
			{
				peak = yval;
				peakLoc = curveData.centre;
			}
debugP("y peak = "+pNum(yval, precision)+" at theme "+(curveData.centre+1));
}
/*
else	// we can't wwork backwards lacking a trough max value!
{	// trough
			if ( yval > trough )
			{
				trough = yval;
				troughLoc = curveData.centre;
			}

}
*/
		}	// done all functions in record
		j += numFunctions+3;
		i++;
		if ( peakLoc == -1 )
debugP("couldn't find peak");
		else
{
debugP("peak max at theme "+(peakLoc+1));
			numPrimaryPeaks[peakLoc] += 1;
}
/*
		if ( troughLoc == -1 )
debugP("couldn't find peak");
		else
{
debugP("trough at theme "+(troughLoc+1));
			numPrimaryTroughs[troughLoc] += 1;
}
*/
		}	// done all records
debugP("theme primary peaks");
for (n=0; n < 12; n++)
debugP((n+1)+" "+numPrimaryPeaks[n]);
/*
var stats = { mean: 0, sd: 0, sError: 0, skew: 0 };
statAnalyse ( numPrimaryPeaks, 0, 12, stats, 0 );
debugP("primary peak stats");
debugP("mean "+pNum(stats.mean, precision)+" std. dev. "+pNum(stats.sd, precision));
/*
/*
stats.mean = 0;
stats.sd = 0;
stats.sError = 0;
stats.skew = 0;
debugP("theme primary troughs");
for (n=0; n < 12; n++)
debugP((n+1)+" "+numPrimaryTroughs[n]);
statAnalyse ( numPrimaryTroughs, 0, 12, stats, 0 );
debugP("primary trough stats");
debugP("mean "+pNum(stats.mean, precision)+" std. dev. "+pNum(stats.sd, precision));
*/
}

function ordinate ( xLoc, curveData, dataSize )
{
	var ixL, ixR, x, y, y1, y2;
	x = xLoc;
	ixL = elementInBounds ( curveData.centre-curveData.distL, dataSize);
	ixR = elementInBounds ( curveData.centre + curveData.distR, dataSize);
	 if ( x > curveData.centre )
			x -= curveData.centre;
	 else x -= ixL;

	 if ( x < 0 )
		 x += 12;

	y = 0;
	y0 = curveData.aQ * Math.pow ( x, 2 ) + curveData.bQ * x + curveData.cQ;	// quadratic
	
	if ( x < curveData.centre )
	{
		y1 = curveData.m0 * x + curveData.c0;	// LHS fn
		y2 = 0;
	}
	else
	{
		y2 = curveData.m1 * ( x - curveData.distL ) + curveData.c1;	// RHS fn
		y1 = 0;
	}
var yAvg;
	y = y0 * ( 1 - fnRatio ) + ( y1 + y2 ) * fnRatio;	// fnRatio is a global in engine.js!
yAvg = y;
	y = y0 * ( 1 - (1/fnRatio) ) + ( y1 + y2 ) * (1/fnRatio);
yAvg = (y + yAvg) * 0.5;
return yAvg;
} 

			

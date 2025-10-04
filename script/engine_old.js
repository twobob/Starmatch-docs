/* ********************************************************************************* */
/* ****************************** IN MEMORIAM ************************************** */
/* **************************** Claudius Ptolemy *********************************** */
/* ********************************************************************************* */

/* ********************************************************************************** */

 /*   This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

var DEBUG = 0;	// allow function expression printout
var rp;	// ruling planet ( number )
var mcSignRuler;

var aspects = [0,180,150,120,90,60,45,30];	// FIXED
var totalAspects;
var af = [1,0.5,0.3727,0.3333,0.3536,0.2357,0.25,0.1667];	// FIXED
/* aoIndex = 0 for default (offset into array): USER-EXTENSIBLE */
/* Note: aspect orb sets 3,4,5 from 'Astrotheme.com: natal, synastry, transit */
/* N.b. orbs traditionally have factors, applying to the planets not the aspects */
/* orbType: 0 = aspect orbs (modern), 1 = planet orbs (traditional) */
var orbType = 0;
/* planet orbs: Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto */
/* the last three are modern additions */
var po = [[15,12,7,7,7,9,9,5,5,5],[17,12.5,7,8,8,12,10,5,5,5]];	// USER-EXTENSIBLE
/* poIndex = 0, 1 (Lilly, al-Biruni) */
var poIndex = 0;
var ao = [[9,9,3,7,7,5,3,3],[9,9,5,9,9,6,2,3],[10.8,2.2,10.0,8.3,7.5,5.7,2.5,1.5],[10,8,2,6,6,4.5,1,1],[2.6,2.5, 1.5,2.3,2.3,1.3,1,1]];	// USER-EXTENSIBLE
var aoIndex = 0;
/* Traditional factors (ancient/modern): [signNum[ruler],[exaltation],[detriment],[fall]] */
/* tfIndex = 0 for ancient, 4, 8 for modern (offset into array): USER-EXTENSIBLE */
var tf = [[4,3,2,1,0,2,3,4,5,6,6,5],[0,1,-1,5,-1,2,6,-1,-1,4,-1,3],[3,4,5,6,6,5,4,3,2,1,0,2],[6,-1,-1,4,-1,3,0,1,-1,5,-1,2],[4,3,2,1,0,2,3,9,5,6,7,8],[9,1,-1,5,8,2,6,7,-1,4,-1,3],[3,9,5,6,7,8,4,3,2,1,0,2],[6,7,-1,4,-1,3,9,1,-1,5,8,2],[4,3,2,1,0,2,3,9,5,6,7,8],[0,1,-1,8,5,2,6,7,-1,4,-1,9],[3,9,5,6,7,8,4,3,2,1,0,2],[6,7,-1,4,-1,9,0,1,-1,8,-1,2]];	// USER-EXTENSIBLE
var tfIndex = 0;
/* Traditional factors:  [signNum[polarity],[triplicity],[quadruplicity]] */
var ptq = [[1,0,1,0,1,0,1,0,1,0,1,0],[0,1,2,3,0,1,2,3,0,1,2,3],[0,1,2,0,1,2,0,1,2,0,1,2]];	// FIXED
/* Theme values ( t[n] in algorithm ) */
var theme = [0,0,0,0,0,0,0,0,0,0,0,0];	// VARIABLE, RESERVED
/* numAspects	[conjunction, opposition, trine, square, sextile, semi-square, semi-sextile] */
var numAspects = [0,0,0,0,0,0,0,0];	// number of aspect types. VARIABLE, RESERVED
/* numTradFactors [+ve, -ve, fire, earth, air, water, cardinal, fixed, mutable] */
var numTradFactors = [0,0,0,0,0,0,0,0,0];	// totals for traditional factors. VARIABLE, RESERVED
/* tfDominant: [polarity, triplicity, quadruplicity] - dominant pol., trip., quad. or -1 */
var tfDominant = [0,0,0];	// VARIABLE, RESERVED
/* psRT planet strength [planetNum,...]  - arbitrary value for contribution to event occurrence */
/* if psRT[n] set to zero, cancels the relevant planet's effect */
var precessionFlag;	// true, precess position data before theme calculation
var precessionVal = -130;	// Hipparchus
var precessedTheme = [0,0,0,0,0,0,0,0,0,0,0,0];	// VARIABLE, RESERVED
var nativityYear;
var orbValue;
var systemPlanets;	// either 7 or 10, depending on ancient or modern trad factors: VARIABLE, RESERVED
var method;	// default fit type
var chartType = 0;	// 0: per-theme; 1: 3-theme moving avg.
var avgThemeVal = [0,0,0,0,0,0,0,0,0,0,0,0];	// VARIABLE, RESERVED
var totalThemes = 12; // n um themes used in analysis - equal house method

var curveNameSet = [];
var precision = 5; 	// default 5 d.p.
var saQ, sbQ, scQ, sdQ, seQ, taQ, tbQ, tcQ, tdQ, teQ;	// only required by Qdump for testing
var fnRatio = 1/Math.sqrt(2);	// ratio of linear to quadratic average ( 0 = quadratic, 1 = linear );
var sigmaFactor = 1;	// 85%
var ctWeighting = 0.7;	// composite theme weighting ('blur' factor for type1 (xProfile) charts)
var reductionFactorDebilityFall = 1;	// i.e. no reduction
var smamOnly = 0;

var psRT = [1,1,0.6,0.5,0,0,0,0,0,0,1,1];	// 125

var afRescale = 1;	// smam only, 80%

for (n = 0; n < af.length; n++ )
	af[n] *= afRescale;

	function signOf ( value )
	{
		if ( value < 0 )
			return -1;
		if ( value > 0 )
			return 1;
		return 0;
	}

	function minOf ( x, y )
	{
		if ( x < y )
			return x;
		else
			return y;
	}
	
	function maxOf ( x, y )
	{
		if ( x > y )
			return x;
		else
			return y;
	}

	function fitLinearEq ( xArray, yArray, index, numValues )
	{
		var dataL = { slope: 0, intercept: 0 };
		var meanComX = calcMean ( xArray, index, numValues, 1 );
		var meanComY = calcMean ( yArray, index, numValues, 0 );
		dataL.slope = calcSlope ( meanComX, meanComY, xArray, yArray, index, numValues, 1 );
		dataL.intercept = calcIntercept ( meanComX, meanComY, dataL.slope );
		return dataL;
	}	

	function calcMean ( array, index, numValues, dereference )
	{
		var mean = 0;
		var length = array.length;
		var n, m, o;

		for ( n = index; n < index+numValues; n++ )
		{
			m = elementInBounds ( n, length );
			if ( dereference )
				o = n - index;
			else
				o = m;

			mean += array[o];
		}
		mean /= numValues;
		return mean;
	}

	function calcSlope ( xMean, yMean, xArray, yArray, index, numValues, dereference )
	{
		var length = xArray.length;
		var slope = 0;
		var nSum = 0;
		var dSum = 0;
		var n, m, o;

		for ( n = index; n < index+numValues; n++ )
		{
			m = elementInBounds ( n, length );
			if ( dereference )
				o = n - index;
			else
				o = m;

			nSum += ( o - xMean ) * ( yArray[m] - yMean );
			dSum += Math.pow ( ( xArray[o] - xMean ), 2 );
		}
		slope = nSum/dSum;
		return slope;
	}

	function calcIntercept ( xMean, yMean, slope )
	{
		var intercept;
		var intercept = yMean - ( slope * xMean );
		return intercept;
	}

	function fitQuadraticEq ( xArray, yArray, index, numValues )
	{	// assumes x/y arrays are the same size and index/numValues within array bounds
		var dataQ = { aQ: 0, bQ: 0, cQ: 0, dQ: 0, eQ: 0, x0Q: 0, x1Q: 0 };
		var aQ, bQ, cQ, dQ, eQ, x0Q, x1Q;
		var n;
		var Sx = 0;
		var Sy = 0;
		var SxPow2 = 0;
		var SxPow3 = 0;
		var SxPow4 = 0;
		var Sxy = 0;
		var SxPow2y = 0;
		var scale = 1/numValues;

		for ( n = index; n < index+numValues; n++ )
		{
			Sx += xArray[n];
			Sy += yArray[n];
			SxPow2 += Math.pow ( xArray[n], 2 );
			SxPow3 += Math.pow ( xArray[n], 3 );
			SxPow4 += Math.pow ( xArray[n], 4 );
			Sxy += ( xArray[n] * yArray[n] );
			SxPow2y += Math.pow ( xArray[n], 2 ) * yArray[n];
		}
		var sXX;
		var sXY;
		var sXXPow2;
		var sXPow2Y;
		var sXPow2XPow2;
		var d0;
		var x0Q;		// roots, currently unused
		var x1Q;
		sXX = SxPow2 - ( Math.pow ( Sx, 2 ) * scale );
		sXY = Sxy - ( Sx * Sy * scale );
		sXXPow2 = SxPow3 - ( ( SxPow2 * Sx ) * scale );
		sXPow2Y = SxPow2y - ( ( SxPow2 * Sy ) * scale );
		sXPow2XPow2 = SxPow4 - ( Math.pow ( SxPow2, 2 ) * scale );

		aQ = ( ( sXPow2Y * sXX ) - ( sXY * sXXPow2 ) ) / ( ( sXX * sXPow2XPow2 ) - Math.pow ( sXXPow2, 2 ) );
		bQ = ( ( sXY * sXPow2XPow2 ) - ( sXPow2Y * sXXPow2 ) ) / ( ( sXX * sXPow2XPow2 ) - Math.pow ( sXXPow2, 2 ) );
		d0 = bQ * bQ;
		cQ = ( Sy * scale ) - ( bQ * ( Sx * scale ) ) - ( aQ * ( SxPow2 * scale ) );
		x0Q = ( -bQ + Math.sqrt ( d0 - ( 4 * aQ * cQ ) ) ) / ( 2 * aQ );
		x1Q = ( -bQ - Math.sqrt ( d0 - ( 4 * aQ * cQ ) ) ) / ( 2 * aQ );
		eQ = cQ - ( d0 / ( 4 * aQ ) );	// predicted extremum

		var d1 = 4 * aQ * ( cQ - eQ );
		var d2 = Math.sqrt ( Math.abs ( d0 - d1) );
		var d3 = -bQ + d2;
		var d4 = 2 * aQ;
		dQ = d3 / d4;	// (float!) location of extremum
		
		dataQ.aQ = aQ;
		dataQ.bQ = bQ;
		dataQ.cQ = cQ;
		dataQ.dQ = dQ;
		dataQ.eQ = eQ;
		dataQ.x0Q = x0Q;
		dataQ.x1Q = x1Q;
		return dataQ;
	}

	function curveAnalyse ( xArray, yArray, centrePos, curveBounds )
	{	// xArray and yArray must have same size
		var curveData = { aQ: 0, bQ: 0,  cQ: 0, m0: 0, c0: 0, m1: 0, c1: 0, centre: 0, distL:0, distR: 0 };
		var ixL = centrePos - curveBounds.distL
		var ixR = centrePos + curveBounds.distR
		var themesMax = xArray.length;
		var xOffset = elementInBounds ( ixL, themesMax );
		var qData, lData;
		var transform = [];
		curveData.distL = curveBounds.distL;
		curveData.distR = curveBounds.distR;
		transform = transposeArray ( ixL, ixR, yArray );	// so 'feature' start is in array[0]
		qData = fitQuadraticEq ( xArray, transform, 0, transform.length );
if ( DEBUG )
{
debugP("x element of ["+elementInBounds( ixL, themesMax )+", "+elementInBounds ( ixR, themesMax )+"], centre "+centrePos); // _signed_ min, max usable x values in fn

debugP("Q fn: y = "+pNum(qData.aQ, precision)+"((x-"+xOffset+")^2)+"+pNum(qData.bQ, precision)+"(x-"+xOffset+")+"+pNum(qData.cQ, precision));
}
		curveData.aQ = qData.aQ;
		curveData.bQ = qData.bQ;
		curveData.cQ = qData.cQ;
		curveData.xOff = centrePos - curveBounds.distL;
		curveData.centre = centrePos;
		lData = fitLinearEq ( xArray, transform, 0, curveBounds.distL+1 );
if ( DEBUG )
debugP("L fn: y = "+pNum(lData.slope, precision)+"(x-"+xOffset+")+"+pNum(lData.intercept, precision));
		curveData.m0 = lData.slope;
		curveData.c0 = lData.intercept;
		lData = fitLinearEq ( xArray, transform, curveBounds.distL, curveBounds.distR+1 );
if ( DEBUG )
debugP("R fn: y = "+pNum(lData.slope, precision)+"(x-"+centrePos+")+"+pNum(lData.intercept, precision));
		curveData.m1 = lData.slope;
		curveData.c1 = lData.intercept;
		return curveData;
	}

	function predictValues ( curveData, referenceArray)
	{	// index and numValues of predicted data must be within bounds of referenceArray
		var v, x, y, y0, y1, y2;
		var aQ = curveData.aQ;
		var bQ = curveData.bQ;
		var cQ = curveData.cQ;
		var xQ = curveData.xQ;
		var m0 = curveData.m0;
		var c0 = curveData.c0;
		var m1 = curveData.m1;
		var c1 = curveData.c1;
		var xOffset = curveData.centre - curveData.distL;
		var numValues = curveData.distL+curveData.distR;

		if ( numValues < referenceArray.length )
			numValues++;		// include RH endpoint value (default if numValues = refArray.length)
		var total = xOffset+numValues;
		debugP("Loc.___Predict___Real____% error");
		var index = 0;	// calculated fn values always refer start theme x value to zero

		for ( v = xOffset; v < total; v++ )
		{
			x = elementInBounds ( v, referenceArray.length );
			y0 = aQ * Math.pow ( index, 2 ) + bQ * index + cQ;	// quadratic
			if ( v < curveData.centre )
			{
				y1 = m0 * index + c0;	// LHS fn
				y2 = 0;
				}
			else
			{
				y2 = m1 * ( index - curveData.distL ) + c1;
				y1 = 0;
			}
			index++;
			y = y0 * ( 1 - fnRatio ) + ( y1 + y2 ) * fnRatio;
			y0 = referenceArray[x];
			debugP("x = "+x+" y = "+pNum(y, precision)+" "+pNum(referenceArray[x], precision)+" "+pNum(100*(1-(y0/y)), precision));
		}
	}

	function transposeArray ( indexL, indexR, sourceArray )
	{	// indices (signed) are relative to centre (extremum) value
		var n, m;
		var transpose = [];

		if ( ( indexR - indexL ) < sourceArray.length )
			indexR++;	// include both endpoints if not using entire array

		for ( n = indexL; n < indexR; n++ )
		{
			m = elementInBounds ( n, sourceArray.length );
			transpose.push( sourceArray[m] );
 		}
		return transpose;
	}

	function elementInBounds ( index, length )
	{
		var n, m;

		n = index;
		m = n; 	// default

		if ( n < 0 )
			m = length + n;
		else
			if ( n >= length )
				m = n -  length;
		return m;
	}

	function peakLimits ( array, centrePos, arrayTruthTable )
	{	// allocate truth array and curveBounds values
		var curveBounds = { distL: 0, distR: 0 };
		var distL = 0;
		var distR = 0;
		var length = array.length;
		var k;

		for ( k = 0; k < length; k++ )
			arrayTruthTable[k] = 0;

		arrayTruthTable[centrePos] = 1; // at least
		var peakL, peakR;
		var testValue = array[ centrePos ];
		peakL = centrePos-1;

		if ( peakL == -1 )
			peakL = length-1;

		while ( array[ peakL ] < testValue )
		{

			arrayTruthTable[peakL] = 1;
			testValue = array[ peakL ];
			peakL--;
			distL++;
			
			if ( peakL < 0 )
				peakL = length-1;
		}
		testValue = array[ centrePos ];
		peakR = centrePos+1;
		
		if ( peakR == length )
			peakR = 0;

		while ( array[ peakR ] < testValue )
		{
			arrayTruthTable[peakR] = 1;
			testValue = array[ peakR ];
			peakR++;
			distR++;
			
			if ( peakR > length-1 )
				peakR = 0;
		}
		curveBounds.distL = distL;
		curveBounds.distR = distR;
		return curveBounds;	// number of elements L/R of extremum pos'n
	}

	function troughLimits ( array, centrePos, arrayTruthTable )
	{	// allocate truth array and curveBounds values
		var curveBounds = { distL: 0, distR: 0 };
		var distL = 0;
		var distR = 0;
		var length = array.length;
		var k;

		for ( k = 0; k < length; k++ )
			arrayTruthTable[k] = 0;

		arrayTruthTable[centrePos] = 1;
		var troughL, troughR;
		var testValue = array[ centrePos ];
		troughL = centrePos-1;

		if ( troughL == -1 )
			troughL = length-1;

		while ( array[ troughL ] >  testValue )
		{
			arrayTruthTable[troughL] = 1;
			testValue = array[ troughL ];
			troughL--;
			distL++;

			if ( troughL < 0 )
				troughL = length-1;
		}
		testValue = array[ centrePos ];
		troughR = centrePos+1;

		if ( troughR == length )
			troughR = 0;

		while ( array[ troughR ] > testValue )
		{
			arrayTruthTable[troughR] = 1;
			testValue = array[ troughR ];
			troughR++;
			distR++;
			
			if ( troughR > ( length-1 ) )
				troughR = 0;
		}
		curveBounds.distL = distL;
		curveBounds.distR = distR;
		return curveBounds;	// number of elements L/R of extremum pos'n
	}

	function statAnalyse ( array, start, numValues, stats, sType )
	{
/*
debugP("engine statAnalyse array");
for (n=0; n<12; n++)
	debugP(pNum(array[n], precision));
*/
		var n;
		var mean = 0;
		mean = calcMean ( array, start, numValues );

		var sd = 0;
		var sError = 0;	// std error on mean
		var tmp = 0;

		for ( n = start; n < numValues; n++)
		{
			m = elementInBounds ( n, numValues );
			tmp += Math.pow ( (array[m] -  mean), 2 );
		}
		tmp /= (numValues-sType);	// 0 1 std. dev, 0 we have whole data not a sample
		sd = Math.sqrt ( tmp );
		sError = sd / Math.sqrt( numValues );
		
		var m3 = 0;
		var pow3 = 0;
		var skew = 0;
		tmp = 0;
		for ( n = start; n < numValues; n++ )
		{
			m = elementInBounds ( n, numValues );
			tmp += Math.pow((array[m] -  mean), 3);
		}
		m3 = tmp/numValues;
		tmp = Math.pow(sd, 2);
		pow3 = Math.pow ( tmp, 1.5 );	// power (3/2)
		skew = m3/pow3;
		stats.mean = mean;
		stats.sd = sd;
		stats.sError = sError;
		stats.skew = skew;
		return stats;
	}

	function pNum ( longNum, precision )
	{
	
if ( checkPrecision(longNum) < precision  )	
{ return longNum;  }
        var num;
		var num2;
		var factor;
		factor = Math.pow ( 10, precision );
		num = longNum * factor;
		num2 = Math.round ( num + 0.5 ) / Math.pow ( 10, precision );
		return num2;
	
	}

	function checkPrecision(value) {
		if (!isFinite(value)) return 0;
		var e = 1, p = 0;
		while (Math.round(value * e) / e !== value) { e *= 10; p++; }
		return p;
}
	
	function themeFit ( sArrayVal, tArrayVal, sTpSgn, tTpSgn, statsS, statsT )
	{	// s, t theme values, s, t inflection sign at theme value, s,t limit values (mean +/- SD default)
		var fit = -1;	// invalid
		var sMin, sMax, tMin, tMax;	// limit values

		sMin = statsS.mean - (statsS.sd*sigmaFactor);
		sMax = statsS.mean + (statsS.sd*sigmaFactor);
		tMin = statsT.mean - (statsT.sd*sigmaFactor);
		tMax = statsT.mean + (statsT.sd*sigmaFactor);

		if ( sTpSgn == tTpSgn )	// both either peaks or troughs, result is similarity
		{

			if ( ( ( sArrayVal < sMin ) && ( tArrayVal < tMin ) ) || ( ( sArrayVal > sMax ) && ( tArrayVal > tMax ) ) )	// t-t or p-p
				fit = Math.abs ( sArrayVal - tArrayVal );	// return +ve fit
		}
		else
		{	// either t/p or p/t, result is complementarity
			if ( ( sArrayVal < sMin ) && (tArrayVal < tMin ) )
				fit = Math.abs ( sArrayVal - tArrayVal );
			
			if ( ( sArrayVal < sMin ) && (tArrayVal > tMax ) )
				fit = Math.abs ( sArrayVal - ( tArrayVal - tMax ) );
				
			if ( ( sArrayVal > sMax ) && (tArrayVal < tMin ) )
				fit = Math.abs ( ( sArrayVal - sMax ) - tArrayVal );
			
			if ( ( sArrayVal > sMax ) && (tArrayVal > tMax ) )
				fit = Math.abs ( ( sArrayVal - sMax ) - ( tArrayVal - tMax ) );
		}
		return fit;
	}
/*
	function scaleFactor ( sCentre, tCentre )	// old method for theme _features_
	{
		var sf;
		sf = Math.abs(sCentre-tCentre);
	
		if ( sf > 6 )
			sf = 12 - sf;
	
		sf = Math.cos((sf/6)*(3.14159/2));	// miserable approx to PI/2, should find the js one...
		return sf;
	}
*/
function scaleThemeValue ( sFeatureWidth, tFeatureWidth, maxT, numT, fit )	// scale fit across peak or trough
{	// no protection against div by zero - extremely unlikely
var featureRatio = minOf ( sFeatureWidth, tFeatureWidth ) / maxOf ( sFeatureWidth, tFeatureWidth );
var relThemeScale;
if ( numT == maxT )	// should never happen
{
	relThemeScale = 1;
	debugP("Error: Both charts have identical data, setting relative theme scale to 1");
}
else
	relThemeScale = 1 / ( maxT - numT ); //invert to scale against whole chart
var fitScale;
fitScale = relThemeScale * featureRatio;
scaledFit = fit * fitScale;
	return scaledFit;
}

function featureSize ( limits )
{
	return limits.distL+limits.distR+1;
}
	
	
	function Qdump ( aQ, dQ, eQ )
	{
		debugP("TP "+pNum ( sdQ, 0 ));
		debugP("aQ "+pNum(aQ, precision));
		if ( aQ > 0 )
			debugP("concave");
		if ( aQ < 0 )
			debugP("convex");
		if ( aQ == 0)
			debugP("linear");
		if ( ( dQ < 0 ) || ( dQ > 11 ) )	// note s/b exL but this is too limiting for approx. curve
			debugP("WARNING! extremum seQ "+pNum(eQ, precision)+" location failed: dQ "+pNum(dQ, precision));
	}

	function isAspect ( pos1, pos2, aspect, orb )
	{
		var diff = Math.abs(pos1-pos2);
		diff = ( diff > 180 ? 360-diff : diff );
		if ( ( aspect-orb < diff ) && ( aspect+orb > diff ) )
			return 1;
		return 0;
	}
			
	function relativeAspectStrength ( pos1, pos2, aspect, orb, factor )
	{	// -ve return is not-an-aspect
		var strength = 0;
		var pDiff = Math.abs ( pos1- pos2 );
		pDiff = ( pDiff > 180 ? 360 - pDiff : pDiff );
		strength = 1 - ( Math.abs ( pDiff - aspect ) / orb );
		strength =strength*factor;
		return strength;
	}
	
	function signNum ( pos )
	{
		var value = pos/30-0.5;
		value = ( value<0 ? 0 : value );
		value = Math.round ( value );
		value = ( value>=12 ? -1 : value );
		return value;
	}

	function house ( pos, ascendant )
	{
		var value = ( pos - ascendant ) / 30-0.5;
		value = Math.round ( value );
		value = ( value < 0 ? value + 12 : value );
		value = ( value > 12 ? value - 12 : value );
		return value+1;
	}

	function locateThemeMax ( array, index, numValues )
	{
		var themeMaxValue = 0;
		var themeMaxNum = -1;
		var v;

		for ( v = index; v < index+numValues; v++ )
			if ( array[v] > themeMaxValue )
			{
				themeMaxValue = array[v];
				themeMaxNum = v;
			}
		return themeMaxNum;
	}

	function locateThemeMin ( array, index, numValues )
	{
		var themeMinValue = 1;
		var themeMinNum = -1;
		var v;

		for ( v = index; v < index+numValues; v++ )
			if ( array[v] <= themeMinValue )
			{
				themeMinValue = array[v];
				themeMinNum = v;
			}
		return themeMinNum;
	}

	function compositeThemeValues ( array )
	{
		var t, u, v;
		for ( v = 0; v < 12; v++ )
			avgThemeVal[v] = 0;

		for ( v = 0; v < 12; v++ )
		{
			for ( u = -1; u < 2; u++ )
			{
				t = v+u;
				
				if ( t > 11 )
					t -= 12;
				if ( t < 0 )
					t = 12+t;

				if ( Math.abs ( u ) == 1 )
					avgThemeVal[v] += (ctWeighting * array[t]);
				else
					avgThemeVal[v] += array[t];
			}
			avgThemeVal[v] /= (1 + 2*ctWeighting);
		}
	}

	function allocateObjectValueToTheme ( objectPos, strength )
	{
		var posValue = objectPos%30;

		if ( posValue > 15 )
			return strength * ( 1 - ( ( posValue - 15 ) / 30 ) );
		else
			return strength * ( ( posValue / 30 ) + 0.5 );
	}

	function getThemeValues(Sun,Moon,Mercury,Venus,Mars,Jupiter,Saturn,Uranus,Neptune,Pluto,Ascendant,Midheaven)
 	{
/*
for ( n = 0; n < av.length; n++ )
	af[n] = av[n];	// restore defaults

if ( smamOnly == 1 )
afRescale = 45;
for (n = 0; n < av.length; n++ )
	af[n] *= afRescale;
*/
		var tfSet = tfIndex*0.25;
		function numPlanetsInHouse ( houseNum )
		{
			numPlanets = 0;
			for ( n = 0; n < systemPlanets; n++ )
				if ( house ( planet[n], planet[10] ) == houseNum )
					numPlanets++;
			return numPlanets;
		}

		function numStrongPlanetsInHouse ( houseNum )
		{
			numStrong = 0;
			for ( n = 0; n < systemPlanets; n++ )
				if ( house ( planet[n], planet[10] ) == houseNum )	// planet n in house 1
				{
					if ( tf[tfIndex][signNum(planet[n])] == n )	// ruler?
						numStrong++;
					if ( tf[tfIndex+1][signNum(planet[n])] != -1 )
						if ( tf[tfIndex+1][signNum(planet[n])] == n )	// exalted?
							numStrong++;
				}
			return numStrong;
		}
		
		function numPlanetsInSign ( sign )
		{
			numPlanets = 0;
			for ( n = 0; n < systemPlanets; n++ )
				if ( signNum ( planet[n] ) == sign )
					numPlanets++;
			return numPlanets;
		}
		
		function numStrongPlanetsInSign ( sign )
		{
			numStrong = 0;
			for ( n = 0; n < systemPlanets; n++ )	// for all planets
			{
				if ( signNum ( planet[n] ) == sign )
				{
					if ( tf[tfIndex][sign] == n )	// ruler
						numStrong++;
					if ( tf[tfIndex+1][signNum(planet[n])] != -1 )
						if ( tf[tfIndex+1][sign] == n )	// exalted
							numStrong++;
				}
			}
			return numStrong;
		}
	
		function isMutualReception ( Px )
		{
			var signY = signNum(planet[Px]);
			if ( tf[tfIndex][signY] != Px )
			{
				var signX;
				for ( m = 0; m < systemPlanets; m++ )
				{
					if ( tf[tfIndex][m] == Px )
					{
						signX = m;
						if ( signNum(planet[tf[tfIndex][signY]]) == signX )
							return signNum(planet[tf[tfIndex][signY]]);
					}
				}
			}
			return -1;
		}
	
		/* Precession of equinoxes - shift of 1st. point of Aries (Ras Hammel still corresponds
		in the Hindu system of sidereal astrology) is now in Pisces.
		Does this make an astrological Aries (outside the sidereal system) a Pisces?
		Precession of equinoxes gives great circle of approx. 25772 years. First point of
		Aries defined in 130 BCE by Hipparchus.
		Current first point of Aries is thus 360*(currentYear+130)/25772, or about 0 Pisces.
		*/
		function precession ( year )
		{
			var elapsedYears;
			elapsedYears = Math.abs ( year - precessionVal );
			return elapsedYears;			
		}
		
		function precessPositions ( nativityYear )
		{
			var degrees = precession (nativityYear );
			for ( n = 0; n < 12; n++ )
			{
				planet[n] -= degrees;
				planet[n] = ( planet[n] < 0 ? 360+planet[n] : planet[n] );
				planet[n] = ( planet[n] > 360 ? planet[n]-360 : planet[n] );
			}
		}

		function dualRuler ( rulingPlanet )	// check for doubled sign rulers
		{
			var q;
			var dual = 0;
			for ( q = 0; q < 12; q++ )	// 12 signs in system
				if ( tf[tfIndex][q] == rulingPlanet )
					dual++;

			if ( dual > 1 )
				return 1;
			else
				return 0;
		}	

		function calculateThemeValue ( themeNum, signRuler, rulerWeighting )	// 1 - 12
		{
var weighting = 1;
if  ( dualRuler ( signRuler ) )
	rulerWeighting *= 0.5;
			var themeValue, inMR;
			themeValue = 0;
	
			if ( signRuler != -1 )	// check for sign ruler in House themeNum
			{	// is ruler in House?
				if ( house ( planet[signRuler], planet[10] ) == themeNum )
					themeValue += allocateObjectValueToTheme ( planet[signRuler], ps[signRuler] * rulerWeighting );// ignore spurious weighting
/*
				inMR = isMutualReception ( signRuler );

				if ( inMR != -1 )
				{	// effective conjunction affects theme of both planets involved
if ( ( themeNum - 1 ) != inMR )
{
var inMRweighting =ps[tf[0][inMR]];
themeValue += inMRweighting;
theme[inMR] += inMRweighting;
}
				}
*/
			}

				if ( house ( planet[0], planet[10] ) == themeNum )
					themeValue += allocateObjectValueToTheme ( planet[0], ps[0] );// ignore spurious weighting
				if ( house ( planet[1], planet[10] ) == themeNum )
				themeValue += allocateObjectValueToTheme ( planet[1], ps[1] );
{
				if ( house ( planet[rp], planet[10] ) == themeNum )
					themeValue += allocateObjectValueToTheme ( planet[rp], ps[rp] );
}
if (mcSignRuler != 3)
{
				if ( house ( planet[mcSignRuler], planet[10] ) == themeNum )
					themeValue += allocateObjectValueToTheme ( planet[mcSignRuler], ps[mcSignRuler] );

}

/*
			themeValue += ( numStrongPlanetsInHouse ( themeNum ) > 1 ? weighting : 0 );	// 1 point for each?
			themeValue += ( numPlanetsInHouse ( themeNum ) > 1 ? weighting :  0 );
*/
			if ( signRuler != -1 )	// check for signRuler in sign
				if ( signNum ( planet[signRuler] ) == themeNum-1 )
					themeValue += allocateObjectValueToTheme ( planet[signRuler], ps[signRuler]*rulerWeighting );
				if ( signNum ( planet[0] ) == themeNum-1 )
					themeValue += allocateObjectValueToTheme ( planet[0], ps[0] );
				if ( signNum ( planet[1] ) == themeNum-1 )
					themeValue += allocateObjectValueToTheme ( planet[1], ps[1] );
if ( signNum ( planet[10] ) == themeNum-1 )
	themeValue += allocateObjectValueToTheme ( planet[10], ps[10] );
/*
			themeValue += ( signNum ( planet[10] ) == themeNum-1 ? ps[10]*weighting : 0 );
*/
if ( signNum ( planet[11] ) == themeNum-1 )
{
	themeValue += allocateObjectValueToTheme ( planet[11], ps[11] );
}
/*
			themeValue += ( numStrongPlanetsInSign ( themeNum-1 ) > 0 ? weighting : 0 );
			themeValue += ( numPlanetsInSign ( themeNum-1 ) > 1 ? weighting : 0 );
*/
			for ( n = 0; n < totalAspects; n++ )	// aspect list
			{
				if ( themeNum != 5 )
				{
					if ( orbType == 0 )	// aspect orbs
						orbValue = ao[aoIndex][n];
					else
						orbValue = 0.5*(po[poIndex][0]+po[poIndex][signRuler]);	// half sum of planet orbs
					if ( isAspect ( planet[0], planet[signRuler], aspects[n], orbValue ) )	// signRuler/Sun aspect
themeValue += Math.sqrt((ps[0]^2)+(ps[signRuler]^2))*rulerWeighting* relativeAspectStrength ( planet[0], planet[signRuler], aspects[n], orbValue, af[n] );
				}

				if ( themeNum != 4 )
				{
					if ( orbType == 0 )	// aspect orbs
						orbValue = ao[aoIndex][n];
					else
						orbValue = 0.5*(po[poIndex][1]+po[poIndex][signRuler]);
					if ( isAspect ( planet[1], planet[signRuler], aspects[n], orbValue ) )	// signRuler/Moon aspect
themeValue += Math.sqrt((ps[1]^2)+(ps[signRuler]^2))*rulerWeighting* relativeAspectStrength ( planet[1], planet[signRuler], aspects[n], orbValue, af[n] );
				}

				if ( orbType == 0 )
					orbValue = ao[aoIndex][n];
				else
					orbValue = po[poIndex][signRuler];	// we don't have a planet orb for asc. or MC
				if ( isAspect ( planet[10], planet[signRuler], aspects[n], orbValue ) )	// signRuler/Ascendant aspect
themeValue += Math.sqrt((ps[10]^2)+(ps[signRuler]^2))*rulerWeighting* relativeAspectStrength ( planet[10], planet[signRuler], aspects[n], orbValue, af[n] );
 					if ( isAspect ( planet[11], planet[signRuler], aspects[n], orbValue ) )	// signRuler/Midheaven aspect
themeValue += Math.sqrt((ps[11]^2)+(ps[signRuler]^2))*rulerWeighting* relativeAspectStrength ( planet[11], planet[signRuler], aspects[n], orbValue, af[n] );
			}
			theme[themeNum-1] += themeValue;	// include with alloc'd value
		}

		var fSun = TidyUpAndFloat(Sun);
		var fMoon = TidyUpAndFloat(Moon);
		var fMercury = TidyUpAndFloat(Mercury);
		var fVenus = TidyUpAndFloat(Venus);
		var fMars = TidyUpAndFloat(Mars);
		var fJupiter = TidyUpAndFloat(Jupiter);
		var fSaturn = TidyUpAndFloat(Saturn);
		var fUranus = TidyUpAndFloat(Uranus);
		var fNeptune = TidyUpAndFloat(Neptune);
		var fPluto = TidyUpAndFloat(Pluto);
		var fAscendant = TidyUpAndFloat(Ascendant);
		var fMidheaven = TidyUpAndFloat(Midheaven);
		var numPlanets = 0;
		var numStrong = 0;
		/* Planetary positions */
		var planet = [fSun,fMoon,fMercury,fVenus,fMars,fJupiter,fSaturn,fUranus,fNeptune,fPluto,fAscendant,fMidheaven];

		var m,n,o;	// loop vars.
		var k,tmp;
		var ps = [0,0,0,0,0,0,0,0,0,0,0,0];
		totalAspects = aspects.length;


		for ( n = 0; n < 12; n++ )
		{
			theme[n] = 0;
			avgThemeVal[n] = 0;
			ps[n] = psRT[n];		// reset to initialised values (1)
		}


		systemPlanets = 10;
		if ( tfIndex == 0 )
		{
			systemPlanets = 7;
			for ( n = systemPlanets; n < 10; n++ )
				ps[n] = 0;	// No Uranus, Neptune, Pluto - weighting zero
		}
	
		if ( precessionFlag != 0 )
			precessPositions ( nativityYear  );
	
		for ( n = 0; n< 9; n++ )
			numTradFactors[n] = 0;
		for ( n = 0; n< 3; n++ )
			tfDominant[n] = 0;
		for ( n = 0; n < 8; n++ )
			numAspects[n] = 0;
	
		for ( n = 0; n < systemPlanets; n++ )	// for all planets ancient or modern
		{
			k = signNum ( planet[n] )

			if ( ptq[0][k] == 1 )	// can only be [0,1]
				numTradFactors[0]++;	// +ve sign
			else
				numTradFactors[1]++;	// -ve sign
	
			if ( ptq[1][k] == 0 )
				numTradFactors[2]++;	// fire
			if ( ptq[1][k] == 1 )
				numTradFactors[3]++;	// earth
			if ( ptq[1][k] == 2 )
				numTradFactors[4]++;	// air
			if ( ptq[1][k] == 3 )
				numTradFactors[5]++;	// water
			
			if ( ptq[2][k] == 0 )
				numTradFactors[6]++;	// cardinal
			if ( ptq[2][k] == 1 )
				numTradFactors[7]++;	// fixed
			if ( ptq[2][k] == 2 )
				numTradFactors[8]++;	// mutable
		}

		for ( n = 10; n < 12; n++ )	// Asc., M.C.
		{
			k = signNum ( planet[n] )
	
			if ( ptq[0][k] == 1 )	// can only be [0,1]
				numTradFactors[0]++;	// +ve sign
			else
				numTradFactors[1]++;	// -ve sign
	
			if ( ptq[1][k] == 0 )
				numTradFactors[2]++;	// fire
			if ( ptq[1][k] == 1 )
				numTradFactors[3]++;	// earth
			if ( ptq[1][k] == 2 )
				numTradFactors[4]++;	// air
			if ( ptq[1][k] == 3 )
				numTradFactors[5]++;	// water
			
			if ( ptq[2][k] == 0 )
				numTradFactors[6]++;	// cardinal
			if ( ptq[2][k] == 1 )
				numTradFactors[7]++;	// fixed
			if ( ptq[2][k] == 2 )
				numTradFactors[8]++;	// mutable
		}

		if ( numTradFactors[0] > numTradFactors[1] )	// polarity
			tfDominant[0] = 1;		// +ve dominant
		else
			tfDominant[0] = -1;		// -ve dominant
		tfDominant[1] = 2;		// default fire

		for ( n = 3; n < 6; n++ )
			if ( numTradFactors[n] > numTradFactors[tfDominant[1]] )
				tfDominant[1] = n;

		for ( n = 2; n < 6; n++ )
			if ( n != tfDominant[1] )
				if ( numTradFactors[n] == numTradFactors[tfDominant[1]] )	// no dominant trip.
					tfDominant[1] = -1;

		var tripDominantFactor;
		if ( tfDominant[1] != -1 )
			tripDominantFactor = 1 - ( 1 / ( numTradFactors[tfDominant[1]] ) );	// avoid zero weighting!
		else
			tripDominantFactor = 0;

		tfDominant[2] = 6;		// default cardinal
		for ( n = 7; n < 9; n++ )
			if ( numTradFactors[n] > numTradFactors[tfDominant[2]] )
				tfDominant[2] = n;
	
		for ( n = 6; n < 9; n++ )
			if ( n != tfDominant[2] )
				if ( numTradFactors[n] == numTradFactors[tfDominant[2]] )	// no dominant trip.
					tfDominant[2] = -1;
		var quadDominantFactor;
		if ( tfDominant[2] != -1 )
			quadDominantFactor = 1 - ( 1 / ( numTradFactors[tfDominant[2]] ) );	
		else
			quadDominantFactor = 0;

		if ( orbType == 0 )
		{
			for ( n = 0; n < 10; n++ )	// avoid mutual Asc. / M.C. aspects
				for ( m = n+1; m < 12; m++ )
					if ( n != m )
					{
						if ( !( ( ( systemPlanets == 7 ) && ( n > 6 ) ) ) )
						{
							for ( o = 0; o < aspects.length; o++ )	// aspect list
								if ( isAspect ( planet[n], planet[m], aspects[o], ao[aoIndex][o] ) )
									numAspects[o]++;
						}
					}
		}
		else
		{
			for ( n = 0; n < systemPlanets; n++ )	// aspects between planets
				for ( m = n+1; m < systemPlanets; m++ )
					if ( n != m )
					{
						orbValue = 0.5*(po[poIndex][n]+po[poIndex][m]);
						for ( o = 0; o < aspects.length; o++ )
							if ( isAspect ( planet[n], planet[m], aspects[o], orbValue ) )
								numAspects[o]++;
					}

			for ( n = 0; n < systemPlanets; n++ )	// aspects from planets to Asc. and M.C.
				for ( m = 10; m < 12; m++ )
				{
					orbValue = po[poIndex][n];
					for ( o = 0; o < aspects.length; o++ )
					{
						if ( isAspect ( planet[n], planet[m], aspects[o], orbValue ) )
							numAspects[o]++;
					}
				}
		}
		nAspects = -1;
		var dominantAspect = 0;
		for ( n = 0; n < aspects.length; n++ )	// n is aspect number
		{
			if ( numAspects[n] > nAspects )
			{
				nAspects = numAspects[n];
				dominantAspect = n;
			}
		}

		for ( n = dominantAspect+1; n < aspects.length; n++ )
			if ( numAspects[n] == nAspects )	// found duplicate
			dominantAspect = -1;	// reset dominant aspect location to invalid
		if ( dominantAspect != -1 )
		{
			aspectFactor = 1 - ( 1 / ( nAspects + 1 ) );	// avoid zero weighting!
			aspectFactor *= af[dominantAspect];
		}
		else
			aspectFactor = 0;

		var sign;
		for ( n = 0 ; n < systemPlanets; n++ )	// not Ascendant or Midheaven
		{
			sign = signNum ( planet[n] );
			if ( ( tf[tfIndex][sign] == n ) || ( tf[tfIndex][sign] == n ) )
{
				ps[n] *= reductionFactorDebilityFall;	// reduce contribution by arbitrary factor
}
		}

		rp = tf[tfIndex][signNum ( planet[10] )];
mcSignRuler = tf[tfIndex][signNum(planet[11])];

tripDominantFactor = 0;
quadDominantFactor = 0;
aspectFactor = 0;

		calculateThemeValue ( 1, tf[tfIndex][0], 1 );
		if ( tfDominant[1] == 2 )		// fire dominant? (index to numTradFactors[tfDominant[1]])
			theme[0] += tripDominantFactor;
		if ( tfDominant[2] == 6 )		// cardinal dominant?
			theme[0] += quadDominantFactor;
		if ( dominantAspect == 0 )		// conjunctions dominant
theme[0] += aspectFactor;
		calculateThemeValue ( 2, tf[tfIndex][1], 1 );
		if ( tfDominant[1] == 3 )		// earth dominant?
			theme[1] += tripDominantFactor;
		if ( tfDominant[2] == 7 )		// fixed dominant?
			theme[1] += quadDominantFactor;
		if ( dominantAspect == 7 )		// semi-sextiles dominant
theme[1] += aspectFactor;
		calculateThemeValue ( 3, tf[tfIndex][2], 1 );
		if ( tfDominant[1] == 4 )		// air dominant?
			theme[2] += tripDominantFactor;
		if ( tfDominant[2] == 8 )		// mutable dominant?
			theme[2] += quadDominantFactor;

		if ( dominantAspect == 5 )		// sextiles dominant
theme[2] += aspectFactor;
		calculateThemeValue ( 4, tf[tfIndex][3], 1 );	//original has factor 2
		if ( tfDominant[1] == 5 )		// water dominant?
			theme[3] += tripDominantFactor;
		if ( tfDominant[2] == 6 )		// cardinal dominant?
			theme[3] += quadDominantFactor;
		if ( dominantAspect == 4 )		// squares dominant
theme[3] += aspectFactor;
		calculateThemeValue ( 5, tf[tfIndex][4], 1 );	// original value 2
		if ( tfDominant[1] == 2 )		// fire dominant?
			theme[4] += tripDominantFactor;
		if ( tfDominant[2] == 7 )		// fixed dominant?
			theme[4] += quadDominantFactor;
		if ( dominantAspect == 3 )		// trines dominant
theme[4] += aspectFactor;
		calculateThemeValue ( 6, tf[tfIndex][5], 1 );
		if ( tfDominant[1] == 3 )		// earth dominant?
			theme[5] += tripDominantFactor;
		if ( tfDominant[2] == 8 )		// mutable dominant?
			theme[5] += quadDominantFactor;
		if ( dominantAspect == 2 || dominantAspect == 6 )		// quincunxes or semisquares dominant?
theme[5] += aspectFactor;
		calculateThemeValue ( 7, tf[tfIndex][6], 1 );
		if ( tfDominant[1] == 4 )		// air dominant?
			theme[6] += tripDominantFactor;
		if ( tfDominant[2] == 6 )		// cardinal dominant?
			theme[6] += quadDominantFactor;
		if ( dominantAspect == 1 )		// oppositions dominant
theme[6] += aspectFactor;
		var ruler = tf[tfIndex][7];
/*
		if ( ruler == 9 )	// Pluto, modern ruler of Scorpio
		{
			calculateThemeValue ( 8, ruler, 0.5 );
			calculateThemeValue ( 8, 4, 0.5 );	// add contribution from ancient ruler Mars
			for ( n = 0; n < 7; n++ )
			{
				if ( orbType == 0 )
					orbValue = ao[aoIndex][n];
				else
					orbValue =  0.5*(po[poIndex][4]+po[poIndex][ruler]);

				if ( isAspect ( planet[4], planet[ruler], aspects[n], orbValue ))	// signRuler/Mars aspect
					theme[7] += ps[4]*ps[ruler]*relativeAspectStrength ( planet[4], planet[ruler], aspects[n], orbValue, af[n] );
			}

		}
		else
*/
			calculateThemeValue ( 8, ruler, 1 );	// ruler dependent on tfIndex
		if ( tfDominant[1] == 5 )		// water dominant?
			theme[7] += tripDominantFactor;
		if ( tfDominant[2] == 7 )		// fixed dominant?
			theme[7] += quadDominantFactor;
		if ( dominantAspect == 2 || dominantAspect == 6 )		// quincunxes or semisquares dominant?
theme[7] += aspectFactor;
		calculateThemeValue ( 9, tf[tfIndex][8], 1 );
		if ( tfDominant[1] == 2 )		// fire dominant? (index to numTradFactors[tfDominant[1]])
			theme[8] += tripDominantFactor;
		if ( tfDominant[2] == 8 )		// mutable dominant?
			theme[8] += quadDominantFactor;
		if ( dominantAspect == 3 )		// trines dominant
		theme[8] += aspectFactor;
		calculateThemeValue ( 10, tf[tfIndex][9], 1 );
		if ( tfDominant[1] == 3 )		// earth dominant?
			theme[9] += tripDominantFactor;
		if ( tfDominant[2] == 6 )		// cardinal dominant?
			theme[9] += quadDominantFactor;
		if ( dominantAspect == 4 )		// squares dominant
theme[9] += aspectFactor;
		var ruler = tf[tfIndex][10];
/*
		if ( ruler == 7 )	// Uranus, modern ruler of Aquarius
		{
			calculateThemeValue ( 11, ruler, 0.5 );
			calculateThemeValue ( 11, 6, 0.5 );	// add contribution from ancient ruler Saturn
			for ( n = 0; n < 7; n++ )
			{
				if ( orbType == 0 )
					orbValue = ao[aoIndex][n];
				else
					orbValue =  0.5*(po[poIndex][4]+po[poIndex][ruler]);

				if ( isAspect ( planet[6], planet[ruler], aspects[n], orbValue ))	// signRuler/Mars aspect
					theme[10] += ps[6]*ps[ruler]*relativeAspectStrength ( planet[6], planet[ruler], aspects[n], orbValue, af[n] );
			}
		}
		else
*/
			calculateThemeValue ( 11, ruler, 1 );	// just use ancient ruler Saturn
		if ( tfDominant[1] == 4 )		// air dominant?
			theme[10] += tripDominantFactor;
		if ( tfDominant[2] == 7 )		// fixed dominant?
			theme[10] += quadDominantFactor;
		if ( dominantAspect == 5 )		// sextiles dominant
theme[10] += aspectFactor;
		var ruler = tf[tfIndex][11];
/*
		if ( ruler == 8 )	// Neptune, modern ruler of Pisces
		{
			calculateThemeValue ( 12, ruler, 0.5 );
			calculateThemeValue ( 12, 5, 0.5 );	// add contribution from ancient ruler Jupiter
			for ( n = 0; n < 7; n++ )
			{
				if ( orbType == 0 )
					orbValue = ao[aoIndex][n];
				else
					orbValue =  0.5*(po[poIndex][4]+po[poIndex][ruler]);

				if ( isAspect ( planet[5], planet[ruler], aspects[n], orbValue ))	// signRuler/Mars aspect
					theme[11] += ps[5]*ps[ruler]*relativeAspectStrength ( planet[5], planet[ruler], aspects[n], orbValue, af[n] );
			}

		}
		else
*/
			calculateThemeValue ( 12, ruler, 1 );
		if ( tfDominant[1] == 5 )		// water dominant?
			theme[11] += tripDominantFactor;
		if ( tfDominant[2] == 8 )		// mutable dominant?
			theme[11] += quadDominantFactor;
		if ( dominantAspect == 7 )		// semi-sextiles dominant?
theme[11] += aspectFactor;
var dominantSignValue = 0;
if ( tfDominant[0] == 1 )	// +ve signs dominant
	dominantSignValue = 1 - ( 1 / numTradFactors[0] );
else
	if ( tfDominant[0] == -1 )	// -ve signs dominant
	dominantSignValue = 1 - ( 1 / numTradFactors[1] );

dominantSignValue = 0;
		if ( tfDominant[0] == 1 )	// + signs dominant
		{
			n = 0;
			while ( n < 12 )
			{
theme[n] += dominantSignValue/6;
				n += 2;
			}
		}
		else
		{
			if ( tfDominant[0] == -1 )	// - signs dominant
			{
			  n = 1;
			  while ( n < 12 )
			  {
theme[n] += dominantSignValue/6;
				 n += 2;
			  }
			}
		}
		
		if ( precessionFlag  != 0 )
		{	// the First Point of Aries is not well-defined
			m = Math.round(precession ( nativityYear  ) / 360  + 0.5 ); 	// also flaky
			for ( n = 0; n < 12; n++ )
			{
				k = ( n-m < 0 ? 12-m : n-m );
				precessedTheme[k] = theme[n];
			}
		}
	}

	method = 1;	// specify default profile method

	function xProfile ( sArrayT, tArrayT )
	{
/*
debugP("engine.js: xProfile method "+method);
debugP("engine: sArray");
for ( n= 0 ; n < 12; n++ )
	debugP(pNum(sArrayT[1][n], precision));
debugP("engine: tArray");
for ( n= 0 ; n < 12; n++ )
	debugP(pNum(tArrayT[1][n], precision));
*/
		if (method>=3)
			return;	// trap erroneous call if made

		var i, j, k;
		var n;
		var pCoincidence;
		var sArrayType0 = [[0,1,2,3,4,5,6,7,8,9,10,11],[0,0,0,0,0,0,0,0,0,0,0,0]];
		var tArrayType0 = [[0,1,2,3,4,5,6,7,8,9,10,11],[0,0,0,0,0,0,0,0,0,0,0,0]];
		var sArrayType1 = [[0,1,2,3,4,5,6,7,8,9,10,11],[0,0,0,0,0,0,0,0,0,0,0,0]];
		var tArrayType1 = [[0,1,2,3,4,5,6,7,8,9,10,11],[0,0,0,0,0,0,0,0,0,0,0,0]];
		var sThemes = [0,0,0,0,0,0,0,0,0,0,0,0];
		var tThemes = [0,0,0,0,0,0,0,0,0,0,0,0];
    curveDataSet = []; // GLOBAL STORE FOR COEFFICIENTS
    curveNameSet = [];
    inflectionPoints = [];

		var arraySize = sThemes.length;	// arrays Themes, Type0/1 must be of same dimension

		for ( n = 0; n < arraySize; n++ )
		{
			sArrayType0[1][n] = sArrayT[1][n];
			tArrayType0[1][n] = TidyUpAndFloat(tArrayT[1][n]);
			sArrayType1[1][n] = sArrayT[1][n];
			tArrayType1[1][n] = TidyUpAndFloat(tArrayT[1][n]);
		}

		var tMinS;
		var tMinT;
		var tMaxS;
		var tMaxT;
		var tPP, tTT;	// diff between peaks/troughs
		compositeThemeValues ( sArrayType1[1] );
		for ( n = 0 ; n < arraySize; n++ )
			sArrayType1[1][n] = avgThemeVal[n];

		compositeThemeValues ( tArrayType1[1] );
		for ( n = 0 ; n < arraySize; n++ )
			tArrayType1[1][n] = avgThemeVal[n];
		var sArrayTp = [[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0]];
		var tArrayTp = [[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0]];
		var sArray = [];
		var tArray = [];
		var data = { fitC: 0, fitS: 0, pCoincidence: 12, sInflections: 0, sArrayTp, sArray, tInflections: 0, tArrayTp, tArray };
		var sStats = { mean: 0, sd: 0, sError: 0, skew: 0 };
		var tStats = { mean: 0, sd: 0, sError: 0, skew: 0 };
/*
debugP("engine: sArray");
for ( n= 0 ; n < 12; n++ )
	debugP(pNum(sArrayType1[1][n], precision));

debugP("engine: tArray");
for ( n= 0 ; n < 12; n++ )
	debugP(pNum(tArrayType1[1][n], precision));
/*
		statAnalyse ( sArrayType1[1], 0, arraySize, sStats, 0 );
		statAnalyse ( tArrayType1[1], 0, arraySize, tStats, 0 );
/*
debugP("engine");
debugP("s mean "+pNum(sStats.mean, 3));
debugP("s sd "+pNum(sStats.sd, 3));
debugP("t mean "+pNum(tStats.mean, 3));
debugP("t sd "+pNum(tStats.sd, 3));
*/
/*
debugP("xProfile statistics:");
debugP("mean: "+pNum(sStats.mean, precision));
debugP("sigma: "+pNum(sStats.sd, precision));
debugP("limits (+/ 1 sigma): "+pNum((sStats.mean-sStats.sd), precision)+", "+pNum((sStats.mean+sStats.sd), precision));
debugP("");
debugP("mean: "+pNum(tStats.mean, precision));
debugP("sigma: "+pNum(tStats.sd, precision));
debugP("limits (+/ 1 sigma): "+pNum((tStats.mean-tStats.sd), precision)+", "+pNum((tStats.mean+tStats.sd), precision));
debugP("");
*/
		if (method == 2)
		{	// use transformed chart limits to find pCoincidence used in method 2, 'type of fit'
			tMinS = locateThemeMin ( sArrayType1[1], 0, arraySize );
			tMinT = locateThemeMin ( tArrayType1[1], 0, arraySize );
			tMaxS = locateThemeMax ( sArrayType1[1], 0, arraySize );
			tMaxT = locateThemeMax ( tArrayType1[1], 0, arraySize );

			pCoincidence = minOf ( Math.abs ( tMinS - tMaxT ), Math.abs ( tMaxS - tMinT ) );
			tPP = Math.abs ( tMaxS - tMaxT );
			tTT = Math.abs ( tMinS - tMinT );	// currently unused

			if ( Math.abs ( tPP ) < pCoincidence )	// i.e. p/p coincidence closer than p/t
				pCoincidence = - ( tPP+1 )	// translate [0, n] -> [-1, -n]

			data.pCoincidence = pCoincidence;
		}

		var inflectionData;
		inflectionData = turningPoints ( sArrayType1 );
		var numInflectionsS = inflectionData.numInflections;
		sArrayTp = inflectionData.tpArray;
		inflectionData = turningPoints ( tArrayType1 );
		var numInflectionsT = inflectionData.numInflections;
		tArrayTp = inflectionData.tpArray;
		var scale;
		var pFit = 0;
		var featureFit = 0;
var totalSignificantThemes = 0;
		var numSignificantThemes;	// num significant themes (used to determine fit)
var numCommonThemes;	// why is this different? other is significvant themes
var curveWidth = { distL: 0, distR: 0 };
var sFeatureSize, tFeatureSize;
		for ( i = 0; i < numInflectionsS; i++ )
		{
			if ( sArrayTp[1][i] == 1 )	// find trough
			{
				curveWidth = troughLimits ( sArrayType1[1], sArrayTp[0][i], sThemes );
sFeatureSize = featureSize ( curveWidth );
				for ( j = 0; j < numInflectionsT; j++ )
				{
numSignificantThemes = 0;
					if ( tArrayTp[1][j] == 1 )	// find trough
					{
featureFit = 0;
						curveWidth = troughLimits ( tArrayType1[1], tArrayTp[0][j], tThemes );
tFeatureSize = featureSize ( curveWidth );
						for ( k = 0; k < arraySize; k++ )
							if ( sThemes[k] & tThemes[k] )
							{
								pFit = themeFit ( sArrayType1[1][k], tArrayType1[1][k],			sArrayTp[1][i], tArrayTp[1][j], sStats, tStats );
								if ( pFit != -1 )	// pfit -ve no fit found
								{
featureFit += pFit
								
									numSignificantThemes++;
								}
							}
					}
					else
					{
						if ( tArrayTp[1][j] == -1 )	// find peak
						{
featureFit = 0;
							curveWidth = peakLimits ( tArrayType1[1], tArrayTp[0][j], tThemes );
tFeatureSize = featureSize ( curveWidth );
							for ( k = 0; k < arraySize; k++ )
								if ( sThemes[k] & tThemes[k] )
								{
									pFit = themeFit ( sArrayType1[1][k], tArrayType1[1][k], sArrayTp[1][i], tArrayTp[1][j], sStats, tStats );
									if ( pFit != -1 )
									{
featureFit -= pFit;
										numSignificantThemes++;
									}
								}
						}
					}
var numCommonThemes;
if ( numSignificantThemes > 0 )
{
numCommonThemes = 0;
totalSignificantThemes += numSignificantThemes;
for (n = 0; n < 12; n++)
if ( sThemes[n] & tThemes[n] )
	numCommonThemes++;
if ( featureFit < 0 )
	data.fitC += scaleThemeValue ( sFeatureSize, tFeatureSize, totalThemes, numCommonThemes, featureFit );
else
	data.fitS += scaleThemeValue ( sFeatureSize, tFeatureSize, totalThemes, numCommonThemes, featureFit );
}
				}	// end j inflections
			}	// end s trough
			else
			{
				if ( sArrayTp[1][i] == -1 )	// find peak
				{
					curveWidth = peakLimits ( sArrayType1[1], sArrayTp[0][i], sThemes );
sFeatureSize = featureSize ( curveWidth );
					for ( j = 0; j < numInflectionsT; j++ )
					{
numSignificantThemes = 0;
						if ( tArrayTp[1][j] == 1 )	// find trough
						{
featureFit = 0;
							curveWidth = troughLimits ( tArrayType1[1], tArrayTp[0][j], tThemes );
tFeatureSize = featureSize ( curveWidth );
							for ( k = 0; k < arraySize; k++ )
								if ( sThemes[k] & tThemes[k] )
								{
									pFit = themeFit ( sArrayType1[1][k], tArrayType1[1][k], sArrayTp[1][i], tArrayTp[1][j], sStats, tStats);
									if ( pFit != -1 )
									{
featureFit = -pFit;
										numSignificantThemes++;
									}
								}
						}
						else
						{
							if ( tArrayTp[1][j] == -1 )	// find peak
							{
featureFit = 0;
								curveWidth = peakLimits ( tArrayType1[1], tArrayTp[0][j], tThemes );
tFeatureSize = featureSize ( curveWidth );
								for ( k = 0; k < arraySize; k++ )
									if ( sThemes[k] & tThemes[k] )
									{
										pFit = themeFit ( sArrayType1[1][k], tArrayType1[1][k], sArrayTp[1][i], tArrayTp[1][j], sStats, tStats);
										if ( pFit != -1 )
										{
featureFit += pFit;
											numSignificantThemes++;
										}
									}
							}
						}
if ( numSignificantThemes > 0 )
{
numCommonThemes = 0;
totalSignificantThemes += numSignificantThemes;
for (n = 0; n < 12; n++)
if ( sThemes[n] & tThemes[n] )
	numCommonThemes++;
if ( featureFit < 0 )
	data.fitC += scaleThemeValue ( sFeatureSize, tFeatureSize, totalThemes, numCommonThemes, featureFit );
else
	data.fitS += scaleThemeValue ( sFeatureSize, tFeatureSize, totalThemes, numCommonThemes, featureFit );

}
					}

				}
			}
		}
		if ( totalSignificantThemes == 0 )	// no fit data found
		{
			data.fitS = -1;	// force invalid fitC and fitS values
			data.fitC = 1;
debugP("engine fit error: fitC "+data.fitC+" fitS "+data.fitS);
		}

		inflectionData = turningPoints ( sArrayType0 );
		data.sInflections = inflectionData.numInflections;
		data.sArrayTp = inflectionData.tpArray;
		data.sArray = sArrayType0[1];
		inflectionData = turningPoints ( tArrayType0 );
		data.tInflections = inflectionData.numInflections;
		data.tArrayTp = inflectionData.tpArray;
		data.tArray = tArrayType0[1];
/*
		data.sInflections = numInflectionsS;
		data.sArrayTp = sArrayTp;
		data.sArray = sArrayType1[1];
		data.tInflections = numInflectionsT;
		data.tArrayTp = tArrayTp;
		data.tArray = tArrayType1[1];
*/

		return data;	// struct {fitC, fitS, pCoincidence, sInflections, sArrayTp, sArray, tInflections, tArrayTp, tArray }
	}	// end xProfile
	
	function turningPoints ( dataArray )	// dataTypeArray (either natal or transformed) // (themeList, themeData)
	{
		var arraySize = dataArray[0].length;	// arrays must be of same dimension
		var tpArray = [[],[]];
		var n;

		for ( n = 0; n < arraySize; n++ )
		{
			tpArray[0][n] = 0;
			tpArray[1][n] = 0;
		}
		var data = { numInflections: 0, tpArray };
		var linearData;		// linear,eq expression coefficients
		var slope = 0;	// subject/target linearData.slope
		var lastSlope = 0;
		var numInflections = 0;
		if ( dataArray[1][arraySize-1] > dataArray[1][0] )
			lastSlope = -1;
		else
			if ( dataArray[1][arraySize-1] < dataArray[1][0] )
				lastSlope = 1;
		
		for ( i = 0; i < arraySize; i++ )
		{
			linearData = fitLinearEq ( dataArray[0], dataArray[1], i, 2 );
			slope = linearData.slope;

			if ( lastSlope != signOf(slope) )
			{
				tpArray[0][numInflections]= i;
				tpArray[1][numInflections] = signOf(slope);
				numInflections++;
			}
	
			lastSlope = signOf(slope);
		}
		if ( dataArray[1][arraySize-1] > dataArray[1][0] )
			slope = -1;
		else
			if ( dataArray[1][arraySize-1] < dataArray[1][0] )
				slope = 1;
	
		if ( ( lastSlope != slope ) && ( slope != 0 ) )
		{
			tpArray[0][numInflections] = arraySize-1;
			tpArray[1][numInflections] = slope;
			numInflections++;
		}
		data.numInflections = numInflections;
		data.tpArray = tpArray;
		return data;	// numInflections, tpArray
	}

	function fnCoeffs ( inflectionPoints )
	{
		var curveDataSet = [];
		curveNameSet = [];	// reset for record.js
		var curveData = { aQ: 0, bQ: 0,  cQ: 0, m0: 0, c0: 0, m1: 0, c1: 0, centre: 0, distL:0, distR: 0 };
		var i, n, m, numInflections, themeNum, tP;
		var limits = { distL: 0, distR: 0 };
		var truthTable;
		var curveData;
		var themeList = [0,1,2,3,4,5,6,7,8,9,10,11];	// dummy data for curve analysis call
		var numItems = inflectionPoints.length;
		curveDataSet.push(numItems/3);	// number of subjects in current selection
		for ( n = 0; n < numItems; n+=3 )	// number of graphs i.e subjects
		{
			m = n;
			numInflections = inflectionPoints[m];
			curveDataSet.push ( m/3 );	// subject id
			curveDataSet.push ( numInflections );	// no. of fns describing this subject's graph
			curveDataSet.push ( inflectionPoints[m+2].length );	// save dataset size
			var themes = inflectionPoints[m+1][0];
			var tp = inflectionPoints[m+1][1];

			for ( i = 0; i < numInflections; i++ )
			{
				themeNum = themes[i];
				truthTable = [0,0,0,0,0,0,0,0,0,0,0,0];
				if ( tp[i] == 1 )
					limits = troughLimits ( inflectionPoints[m+2], themeNum, truthTable );
				else
					limits = peakLimits ( inflectionPoints[m+2], themeNum, truthTable );

				curveData = curveAnalyse ( themeList, inflectionPoints[m+2], themeNum, limits );
				curveData.distL = limits.distL;
				curveData.distR = limits.distR;
				curveData.centre = themeNum;
				curveDataSet.push ( curveData );
if ( DEBUG )
	predictValues ( curveData, inflectionPoints[m+2] )	// check
			}
		}
		return curveDataSet;	// returned data is coeff set for a given record
	}

function normaliseValues ( array, arraySize )
{
		var themeMax = 0;
		for ( n = 0; n < arraySize; n++ )
			if ( array[n] > themeMax )
				themeMax = array[n];

		for ( n = 0; n < arraySize; n++ )
			array[n] /= themeMax;
}
function TidyUpAndFloat(theValue) {
    return parseFloat(theValue);
}

/*
debugP("Natal statistics");
*/
/*
debugP("engine stats");
debugP("mean: s = "+pNum(sStats.mean, precision)+" t = "+pNum(tStats.mean, precision));
debugP("sigma: s = "+pNum(sStats.sd, precision)+" t = "+pNum(tStats.sd, precision));
*/
/*
debugP("std error on mean: s = "+pNum(sStats.sError, precision)+" t = "+pNum(tStats.sError, precision));
*/
/*
debugP("skew: s = "+pNum(sStats.skew, precision)+" t = "+pNum(tStats.skew, precision));
*/

/*
				if ( Math.abs ( sStats.skew ) > Math.abs ( tStats.skew ) )
debugP("S: tail greater or fatter");
				else
debugP("T: tail greater or fatter");

				if ( sStats.skew < 0 )
debugP("S: L tail dominant");
				else
debugP("S: R tail dominant");

				if ( tStats.skew < 0 )
debugP("T: L tail dominant");
				else
debugP("T: R tail dominant");
*/

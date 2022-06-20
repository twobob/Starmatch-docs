/* ********************************************************************************* */
/* *************************** IN MEMORIAM ************************************* */
/* ************************** Claudius Ptolemy ************************************ */
/* ********************************************************************************* */

// In lieu of copyright
// This method of analysis by summing related occurrences of signifiers I found
// in a work by an American woman astrologer whose name, and book title, I do
// not unfortunately recall, about 3 or 4 decades ago.
// Since it is a sensible and logical approach I have followed it here, with some
// modifications.
// Will 18, 2016.
// wj18@talktalk.net, eighteenwill@gmail.com
//
// Credits:
// Program: twobob
// Engine: Will 18
// Engine version: 7.6.3b
/* ********************************************************************************** */

var rp = -1;	// ruling planet ( number )
// AspectValues a: Conjunction, Opposition, Trine, Square, Sextile, Semi-square, Semi-sextile
var a = [0,180,120,90,60,45,30];	// FIXED
// Aspect Factor af: Conjunction, Opposition, Trine, Square, Sextile, Semi-square, Semi-sextile
// values calculated as _reciprocals_ of:
// 1*sqrt(1), 2*sqrt(1), 3*sqrt(1), 2*sqrt(2), 3*sqrt(2), 2*sqrt(4), 3*sqrt(4) (to sqrt(2^n) I suppose...
var af = [1,0.5,0.3333,0.3536,0.2357,0.25,0.1667];	// FIXED
// Aspect Orb ao:  Conjunction, Opposition, Trine, Square, Sextile, Semi-square, Semi-sextile
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
var ao = [[9,9,7,7,5,3,3],[9,9,9,9,6,2,3],[10.8,10.0,8.3,7.5,5.7,2.5,1.5],[10,8,6,6,4.5,1,1],[2.6,2.5,2.3,2.3,1.3,1,1]];
var aoIndex = 0;
/* Traditional factors (ancient/modern): [signNum[ruler],[exaltation],[detriment],[fall]] */
/* tfIndex = 0 for ancient, 4 for modern (offset into array): USER-EXTENSIBLE */
var tf = [[4,3,2,1,0,2,3,4,5,6,6,5],[0,1,-1,5,-1,2,6,-1,-1,4,-1,3],[3,4,5,6,6,5,4,3,2,1,0,2],[6,-1,-1,4,-1,3,0,1,-1,5,-1,2],
[4,3,2,1,0,2,3,9,5,6,7,8],[0,1,-1,5,8,2,6,7,-1,4,-1,3],[3,9,5,6,7,8,4,3,2,1,0,2],[6,7,-1,4,-1,3,9,1,-1,5,8,2]];	// USER-EXTENSIBLE
var tfIndex = 0;
/* Traditional factors:  [signNum[polarity],[triplicity],[quadruplicity]] */
var ptq = [[1,0,1,0,1,0,1,0,1,0,1,0],[0,1,2,3,0,1,2,3,0,1,2,3],[0,1,2,0,1,2,0,1,2,0,1,2]];	// FIXED
/* Theme values ( t[n] in algorithm ) */
var theme = [0,0,0,0,0,0,0,0,0,0,0,0];	// VARIABLE, RESERVED
/* numAspects	[conjunction, opposition, trine, square, sextile, semi-square, semi-sextile] */
var numAspects = [0,0,0,0,0,0,0];	// number of aspect types. VARIABLE, RESERVED
/* numTradFactors [+ve, -ve, fire, earth, air, water, cardinal, fixed, mutable] */
var numTradFactors = [0,0,0,0,0,0,0,0,0];	// totals for traditional factors. VARIABLE, RESERVED
/* tfDominant: [polarity, triplicity, quadruplicity] - dominant pol., trip., quad. or -1 */
var tfDominant = [0,0,0];	// VARIABLE, RESERVED
/* psRT planet strength [planetNum,...]  - arbitrary value for contribution to event occurrence */
/* if psRT[n] set to zero, cancels the relevant planet's effect */
var psRT = [1,1,1,1,1,1,1,1,1,1,1,1];	// VARIABLE, USER-DEFINABLE
var precessionFlag;	// true, precess position data before theme calculation
var precessionVal = -130;	// Hipparchus
var precessedTheme = [0,0,0,0,0,0,0,0,0,0,0,0];	// VARIABLE, RESERVED
var nativityYear;
var orbValue;
var systemPlanets;	// either 7 or 10, depending on ancient or modern trad factors: VARIABLE, RESERVED
var charWeighting = [];	// VARIABLE, RESERVED
var method;	// default fit
var chartType = 0;	// 0: per-theme; 1: 3-theme moving avg.
var avgThemeVal = [0,0,0,0,0,0,0,0,0,0,0,0];	// VARIABLE, RESERVED

// misc global vars
var precision = 5; 	// default 5 d.p.
//var saQ, sbQ, scQ, sdQ, seQ, taQ, tbQ, tcQ, tdQ, teQ;
//var numC;	// ?
//var numS;
var sdFactor = 1;
var sTP;     // num subject inflections (s/b local!)
var tTP;	// num target inflections (s/b local!)
var fnRatio = 1/Math.sqrt(2);	// ratio of linear to quadratic average ( 0 = quadratic, 1 = linear );
// N.b. in tests typically 3 of the 4 cases for S/T P/T prediction have less error than with fnRatio 0.5

// math. functions
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
		var n; var m; var o;

		// if dereference 0 we use the base index
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
		var n; var m; var o;

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
		// xArray is data location, yArray is value
		var dataQ = { aQ: 0, bQ: 0, cQ: 0, dQ: 0, eQ: 0, x0Q: 0, x1Q: 0 };
		var aQ; var bQ; var cQ; var dQ; var eQ; var x0Q; var x1Q;
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
		sXY = Sxy - ( ( Sx * Sy ) * scale );
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
		// curveData: a, b, c coefficients of quadratic eq. y = ax^2+bx+c
		// m, c coefficients of linear eq. y = mx+c for L and R halves of the quadratic
		var curveData = { aQ: 0, bQ: 0,  cQ: 0, m0: 0, c0: 0, m1: 0, c1: 0, centre: 0, ixL:0, ixR: 0 };
		var themesMax = xArray.length;
		var qData; var lData; var transform;

		curveData.ixL = centrePos - curveBounds.distL;
		curveData.ixR = centrePos + curveBounds.distR;
		transform = transposeArray ( curveData.ixL, curveData.ixR, yArray );	// so 'feature' start is in array[0]
		qData = fitQuadraticEq ( xArray, transform, 0, transform.length );
		curveData.aQ = qData.aQ;
		curveData.bQ = qData.bQ;
		curveData.cQ = qData.cQ;
		curveData.xOff = centrePos - curveBounds.distL;
		curveData.centre = centrePos;
		lData = fitLinearEq ( xArray, transform, 0, curveBounds.distL+1 );
		curveData.m0 = lData.slope;
		curveData.c0 = lData.intercept;
		lData = fitLinearEq ( xArray, transform, curveBounds.distL, curveBounds.distR+1 );
		curveData.m1 = lData.slope;
		curveData.c1 = lData.intercept;
		return curveData;
	}

	function predictValues ( limits, curveData, referenceArray)
	{	// index and numValues of predicted data must be within bounds of referenceArray
		var v; var x; var y; var y0; var y1; var y2;
		var aQ = curveData.aQ;
		var bQ = curveData.bQ;
		var cQ = curveData.cQ;
		var xQ = curveData.xQ;
		var m0 = curveData.m0;
		var c0 = curveData.c0;
		var m1 = curveData.m1;
		var c1 = curveData.c1;
		var xOffset = curveData.ixL;
		var numValues = limits.distL+limits.distR;
		var total; var transformedCentre;
		var centrePos = curveData.centre;
		var avgError = 0;

		if ( numValues < referenceArray.length )
			numValues++;		// include RH endpoint value (default if numValues = refArray.length)

		total = xOffset+numValues;
		transformedCentre = centrePos - xOffset;	// keep offset in fn y2 for clarity

		debugP("pV Loc.__Predict__Real____% error");
		for ( v = xOffset; v < total; v++ )
		{
			x = elementInBounds ( v, referenceArray.length );
			y0 = aQ * Math.pow ( ( v - xOffset ), 2 ) + bQ * ( v - xOffset ) + cQ;	// quadratic
			if ( v >= xOffset ) 
			{
				if ( v <= centrePos )
				{
					y1 = m0 * ( v - xOffset ) + c0;	// LHS fn
					y2 = 0;
				}
				else
				{
					y2 = m1 * ( v - xOffset - transformedCentre ) + c1;	// RHS fn
					y1 = 0;
				}
			}
			// (skewed) average better than RMS (quadratic has greater error with larger skew)
			y = y0 * ( 1 -  fnRatio ) + ( y1 + y2 ) * fnRatio;
			y0 = referenceArray[x];
			avgError +=Math.abs(100*(1-(y0/y)));
			avgError /= numValues;
			debugP("x = "+x+" y = "+pNum(y, precision)+" "+pNum(referenceArray[x], precision)+" "+pNum(100*(1-(y0/y)), precision));
		}
		debugP("avg. error "+pNum(avgError, precision));
	}

	function transposeArray ( indexL, indexR, sourceArray )
	{	// indices (signed) are relative to centre (extremum) value
		var n; var m;
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
		var n; var m;

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
		var peakL; var  peakR;
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
		var troughL; var  troughR;
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

	function statAnalyse ( array, start, numValues, stats )
	{
		// standard statistical analysis
		var n;
		var mean = 0;
		mean = calcMean ( array, start, numValues );
		// Math.sqrt(variance) is population std deviation (sigma)
		var sd = 0;
		var tmp = 0;

		for ( n = start; n < numValues; n++)
		{
			m = elementInBounds ( n, numValues );
			tmp += Math.pow ( (array[m] -  mean), 2 );
		}
		tmp /= numValues;
		sd = Math.sqrt ( tmp );

		// skew (method of moments)
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
		stats.skew = skew;
		return stats;
	}

	function pNum ( longNum, precision )
	{
		var num;
		var num2;
		var factor;
		factor = Math.pow ( 10, precision );
		num = longNum * factor;
		num2 = Math.round ( num +0.5 ) / Math.pow ( 10, precision );
		return num2;
	}

// utility functions
	function themeFit ( sArrayVal, tArrayVal, sTpSgn, tTpSgn, statsS, statsT )
	{	// s, t theme values, s, t inflection sign at theme value, s,t limit values (mean +/- SD default)
		var fit = 0;
		var sMin; var sMax; var tMin; var tMax;	// limit values
		var 
		// these limits are FIXED here by use of stats struct
		sMin = statsS.mean - statsS.sd;
		sMax = statsS.mean + statsS.sd;
		tMin = statsT.mean - statsT.sd ;
		tMax = statsT.mean + statsT.sd ;
	
		if ( sTpSgn == tTpSgn )
		{	// signs both the same hence shapes are similar
			if ( sTpSgn == 1 )	// trough, don't need to check both
			{
				if ( ( sArrayVal < sMin ) && ( tArrayVal < tMin ) )	// ONLY if both below limit
					fit = Math.abs( sArrayVal - tArrayVal );
			}
			else
			{	// sgn of BOTH is -ve, check peak
				if ( ( sArrayVal > sMax ) && ( tArrayVal > tMax ) )	// ONLY if both above limit
					fit = Math.abs( sArrayVal - tArrayVal );
			}
		}
		else
		{	// signs differ so shapes are mutually inverse (complementary)
			if ( ( sTpSgn == 1 ) && ( tTpSgn == -1 ) )
			{	// trough/peak so <limit, >limit
				if ( ( sArrayVal < sMin ) && ( tArrayVal > tMax ) )
				{	// find diff as if in same _range_ as similarity
					fit = -Math.abs ( sArrayVal - ( tArrayVal - tMax ) );
				}
			}
			else
			{
				if ( ( sTpSgn == -1 ) && ( tTpSgn == 1 ) )
				{	// peak/trough so >limit, <limit
					if ( ( sArrayVal > sMax ) && ( tArrayVal < tMin ) )
						fit = -Math.abs ( tArrayVal - ( sArrayVal - sMax ) );
				}
			}
		}
		return fit;
	}
	
	function scaleFactor ( sCentre, tCentre )
	{
		var sf;
		sf = Math.abs(sCentre-tCentre);
	
		if ( sf > 6 )
			sf = 12 - sf;
	
		sf = Math.cos((sf/6)*(3.14159/2));	// miserable approx to PI/2, should find the js one...
		return sf;
	}

	function Qdump ( saQ, taQ, sdQ, tdQ, seQ, teQ )
	{
		debugP("s TP "+pNum ( sdQ, 0 )+" t TP "+pNum ( tdQ, 0 ) );
		debugP("saQ "+pNum(saQ, precision)+" taQ "+pNum(taQ, precision));
		if ( saQ > 0 )
			debugP("subject concave");
		if ( saQ < 0 )
			debugP("subject convex");
		if ( saQ == 0)
			debugP("subject linear?");
		if ( taQ > 0 )
			debugP("target concave");
		if ( taQ < 0 )
			debugP("target convex");
		if ( taQ == 0)
			debugP("target linear?");
		if ( ( sdQ < 0 ) || ( sdQ > 11 ) )	// note s/b exL but this is too limiting for approx. curve
			debugP("WARNING! extremum seQ "+pNum(seQ, precision)+" location failed: sdQ "+pNum(sdQ, precision));
		if ( ( tdQ < 0 ) || ( tdQ > 11 ) )	// note s/b exL but this is too limiting for approx. curve
			debugP("WARNING! extremum teQ "+pNum(teQ, precision)+" location failed: tdQ "+pNum(tdQ, precision));
	}

	function isAspect ( pos1, pos2, aspect, orb )
	{
		var diff = Math.abs(pos1-pos2);
		diff = ( diff > 180 ? 360-diff : diff );
		// if diff within aspect-orb and aspect+orb, is aspect
		if ( ( aspect-orb < diff ) && ( aspect+orb > diff ) )
			return 1;
		return 0;
	}
			
	function aspectStrength ( pos1, pos2, aspect, orb, factor )
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

	function relStrength ( sArrayT, tArrayT )
	{
		var v;
		for ( v = 0; v < 12; v++ )
		{
			charWeighting[0] = 0;
			charWeighting[1] = 0;
		}

		for ( v = 0; v < 12; v++ )
		{
			charWeighting[0] += sArrayT[v];	// add theme values
			charWeighting[1] += tArrayT[v];
		}
		// normalise results
		var total  = charWeighting[0]+charWeighting[1];
		charWeighting[0] /= total;
		charWeighting[1] /= total;	// redundant for graphical purposes
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

	function movingAverage ( array )
	{
		var t; var u; var v;
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
				
				avgThemeVal[v] += array[t];
			}
			avgThemeVal[v] /= 3;
		}
	}

	
	
	
	
	
	
			// functions internal to getThemeValues() But stupidly the JS engine wont swallow it
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

		function calculateThemeValue ( themeNum, signRuler, rulerWeighting )	// 1 - 12
		{
			// split allocation between alternate rulers if modern trad. factgors selected
			var weighting = ( rulerWeighting == 0.5 ? 0.5 : 1 );
			// avoid adding contribution from non-ruler associations twice
			var themeValue; var  inMR;
			themeValue = 0;
	
			// i)  Are any of the following in House themeNum? 
			// ruler of sign themeNum-1, Sun, Moon, Ascendant ruler, a strong
			// planet, two or more planets  (allocate one point ( * weighting ) for each).
			if ( signRuler != -1 )	// check for sign ruler in House themeNum
			{	// is ruler in House?
				themeValue += ( house ( planet[signRuler], planet[10] ) == themeNum ? ps[signRuler]*rulerWeighting : 0 );
				// is the ruler in mutual reception?
				inMR = isMutualReception ( signRuler );
				if ( inMR != -1 )
				{	// effective conjunction affects theme of both planets involved
					themeValue += ps[signRuler];	// add a point - note: we do not consider aspect just sign
					theme[inMR] += ps[tf[0][inMR]];
				}
			}
			if ( themeNum != 5 )	// Sun in House themeNum (not Leo)?
				themeValue += ( house ( planet[0], planet[10] ) == themeNum ? ps[0]*weighting : 0 );
			if ( themeNum != 4 )	// // Moon in House themeNum (not Cancer)?
				themeValue += ( house ( planet[1], planet[10] ) == themeNum ? ps[1]*weighting : 0 );
			// Asc. ruler in House themeNum?
			if ( rp != - 1)
				themeValue += ( house ( planet[rp], planet[10] ) == themeNum ? ps[rp]*weighting : 0 );
			// any strong planets? Add extra point
			// it is  possible for more than 1 planet to be in exaltation, depending on
			// rules in tf[planetNum, 1]
			themeValue += ( numStrongPlanetsInHouse ( themeNum ) > 1 ? weighting : 0 );	// 1 point for each?
			// 2 or more planets in house 1?
			themeValue += ( numPlanetsInHouse ( themeNum ) > 1 ? weighting :  0 );

			// ii)  Are any of the following in sign themeNum?  Sign ruler, Sun, Moon, Ascendant
			// a strong planet, two or more planets?
			if ( signRuler != -1 )	// check for signRuler in sign
				themeValue += ( signNum ( planet[signRuler] ) == themeNum-1 ? ps[signRuler]*rulerWeighting : 0 );

			if ( themeNum != 5 )
				themeValue += ( signNum ( planet[0] ) == themeNum-1 ? ps[0]*weighting : 0 );
			if ( themeNum != 4 )
				themeValue += ( signNum ( planet[1] ) == themeNum-1 ? ps[1]*weighting : 0 );
			themeValue += ( signNum ( planet[10] ) == themeNum-1 ? weighting : 0 );
			// strong planets in sign include both rp and exalted
			themeValue += ( numStrongPlanetsInSign ( themeNum-1 ) > 0 ? weighting : 0 );
			// 2 or more planets in sign themeNum-1?
			themeValue += ( numPlanetsInSign ( themeNum-1 ) > 1 ? weighting : 0 );

			// iii) Check for: sign ruler aspecting the Sun, Moon, Ascendant (add
			// strength of the aspect).
			for ( n = 0; n < 7; n++ )	// aspect list
			{
				if ( themeNum != 5 )
				{
					if ( orbType == 0 )	// aspect orbs
						orbValue = ao[aoIndex][n];
					else
						orbValue = 0.5*(po[poIndex][0]+po[poIndex][signRuler]);	// half sum of planet orbs
					if ( isAspect ( planet[0], planet[signRuler], a[n], orbValue ) )	// signRuler/Sun aspect
						themeValue += ps[0]*ps[signRuler]*weighting * rulerWeighting * aspectStrength ( planet[0], planet[signRuler], a[n], orbValue, af[n] );
				}

				if ( themeNum != 4 )
				{
					if ( orbType == 0 )	// aspect orbs
						orbValue = ao[aoIndex][n];
					else
						orbValue = 0.5*(po[poIndex][1]+po[poIndex][signRuler]);
					if ( isAspect ( planet[1], planet[signRuler], a[n], ao[aoIndex][n] ) )	// signRuler/Moon aspect
						themeValue += ps[1]*ps[signRuler]*weighting * rulerWeighting * aspectStrength ( planet[1], planet[signRuler], a[n], orbValue, af[n] );
				}
				
				if ( isAspect ( planet[10], planet[signRuler], a[n], ao[aoIndex][n] ) )	// signRuler/Ascendant aspect
				{
					if ( orbType == 0 )	// aspect orbs
						orbValue = ao[aoIndex][n];
					else
						orbValue = po[poIndex][signRuler];	// we don't have a planet orb for asc.
					themeValue +=ps[signRuler]* weighting * aspectStrength ( planet[10], planet[signRuler], a[n], orbValue, af[n] );
				}

//				if ( themeNum != 1 )	// don't consider MC in Aries?
//				{
					if ( orbType == 0 )	// aspect orbs
						orbValue = ao[aoIndex][n];
					else
						orbValue = po[poIndex][signRuler];	// we don't have a planet orb for MC
					if ( isAspect ( planet[11], planet[signRuler], a[n], ao[aoIndex][n] ) )	// signRuler/Midheaven aspect
						themeValue += ps[signRuler]*weighting * aspectStrength ( planet[11], planet[signRuler], a[n], orbValue, af[n] );
//				}

			}
			theme[themeNum-1] += themeValue;	// allow for multiple calls
		}
	
	
	
	
	
	
	
	
	
	
	
	
	// main algorithm
	function getThemeValues(Sun,Moon,Mercury,Venus,Mars,Jupiter,Saturn,Uranus,Neptune,Pluto,Ascendant,Midheaven)
	{

	
		// Create f prefixed usable values
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
	
		// algorithm main starts here
		var m; var n; var o;	// loop vars.
		var k; var tmp;
		var ps = [0,0,0,0,0,0,0,0,0,0,0,0];
		// initialise theme array

		for ( n = 0; n < 12; n++ )
		{
			theme[n] = 0;
			avgThemeVal[n] = 0;
			ps[n] = psRT[n];		// reset to initialised values
		}

		systemPlanets = 10;
		// no modern planets if factors are traditional not modern (or user defined)
		if ( tfIndex == 0 )
		{
			for ( n= 9; n < 12; n++ )
				ps[n] = 0;	// No Uranus, Neptune, Pluto - weighting zero
			systemPlanets = 7;
		}
	
		if ( precessionFlag != 0 )
			precessPositions ( nativityYear  );
	
		for ( n = 0; n< 9; n++ )
			numTradFactors[n] = 0;
		for ( n = 0; n< 3; n++ )
			tfDominant[n] = 0;
		for ( n = 0; n < 8; n++ )
			numAspects[n] = 0;
	
		// find number of polarities, triplicities, quadruplicities
		//-ve, +ve, fire, earth, air, water, card, fix, mut totals
		for ( n = 0; n <12; n++)	// for all planets, Asc., M.C.	
		{
			k = signNum ( planet[n] );
	
			if ( ptq[0][k] == 1 )	// can only be [0,1]
				numTradFactors[0]++;	// +ve sign
			else
				numTradFactors[1]++;	// -ve sign
	
			if ( ptq[1][k] == 0 )
				numTradFactors[2]++;	// fire
			if ( ptq[1][k] ==1 )
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
			tfDominant[0] = 0;		// -ve dominant
	
		tfDominant[1] = 0;		// default fire
	
		for ( n = 3; n < 6; n++ )
			if ( numTradFactors[n] > tfDominant[1] )
				tfDominant[1] = numTradFactors[n];
	
		for ( n = 3; n < 6; n++ )
			if ( n != tfDominant[1] )
				if ( numTradFactors[n] == tfDominant[1] )	// no dominant trip.
					tfDominant[1] = -1;
	
		tfDominant[2] = 0;		// default cardinal
		for ( n = 6; n < 9; n++ )
			if ( numTradFactors[n] > tfDominant[2] )
				tfDominant[2] = numTradFactors[n];
	
		for ( n = 6; n < 9; n++ )
			if ( n != tfDominant[2] )
				if ( numTradFactors[n] == tfDominant[2] )	// no dominant trip.
					tfDominant[2] = -1;

		// find number of each aspect
		// for each planet, then for each aspect, if aspect, add 1
		if ( orbType == 0 )
		{
			for ( n = 0; n < 12; n++ )
				for ( m = n+1; m < 12; m++ )
					if ( n != m )
					{
						if ( !( ( ( systemPlanets == 7 ) && ( ( n > 6 ) && ( n < 10 ) ) ) ) )
						{
							for ( o = 0; o < 7; o++ )	// aspect list
								if ( isAspect ( planet[n], planet[m], a[o], ao[aoIndex][o] ) )
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
						for ( o = 0; o < 10; o++ )
							if ( isAspect ( planet[n], planet[m], a[o], orbValue ) )
								numAspects[o]++;
					}

			for ( n = 0; n < systemPlanets; n++ )	// aspects from planets to Asc. and M.C.
				for ( m = 10; m < 12; m++ )
				{
					orbValue = po[poIndex][n];
					for ( o = 0; o < 10; o++ )
					{
						if ( isAspect ( planet[n], planet[m], a[o], orbValue ) )
							numAspects[o]++;
					}
				}
		}

		// find dominant aspect (if any)
		tmp = 0;
		for ( n = 0; n < 7; n++ )
			if ( numAspects[n] > tmp )
				tmp = n;
	
		for ( m = 0; m < 7; m++ )
			if ( m != tmp )
			{
				if ( numAspects[m] == tmp )	// no dominant aspect type
					tmp = -1;
			}
		dominantAspect = tmp;
	
		// check for debility/fall in planets and reduce contribution
		for ( n = 0 ; n < systemPlanets; n++ )	// not Ascendant or Midheaven
		{
			var sign = signNum ( planet[n] );
			
			if ( ( tf[3][sign] == n ) || ( tf[4][sign] == n ) )
				ps[n] = 0.5;	// reduce contribution
		}
	
		// Ascendant ruler ( chart ruler )
		rp = tf[tfIndex][signNum ( planet[10] )];
		// ##### end initialisation #####
		
		// Chart analysis - calculate theme values
		// n.b. signs are numbered [0,11], houses are [1,12]
	
		// ##### Theme 1 #####
		calculateThemeValue ( 1, tf[tfIndex][0],1 );
		// Is the chart emphasis on any of the following:  fire, cardinal, conjunctions?
		if ( tfDominant[1] == 0 )		// fire dominant?
			theme[0] += 1;
		if ( tfDominant[2] == 0 )		// cardinal dominant?
			theme[0] += 1;
		if ( dominantAspect == 0 )		// conjunctions dominant
			theme[0] += 1;
		// ##### end Theme 1 #####
	
		// ##### Theme 2 #####
		calculateThemeValue ( 2, tf[tfIndex][1], 1 );
		// Is the chart emphasis on any of the following:  earth, fixed
		if ( tfDominant[1] == 1 )		// earth dominant?
			theme[1] += 1;
		if ( tfDominant[2] == 1 )		// fixed dominant?
			theme[1] += 1;
		theme[1] = ( theme[1] > 3 ? theme[1]+1 : theme[1] );	// this requires max aspect strength!
		// ##### end Theme 2 #####
	
		// ##### Theme 3 #####
		calculateThemeValue ( 3, tf[tfIndex][2], 1 );
		// Is the chart emphasis on any of the following:  air, mutable, sextile
		if ( tfDominant[1] == 2 )		// air dominant?
			theme[2] += 1;
		if ( tfDominant[2] == 2 )		// mutable dominant?
			theme[2] += 1;
		if ( dominantAspect == 4 )		// sextiles dominant
			theme[2] += 1;
		// ##### end Theme 3 #####
	
		// ##### Theme 4 #####
		calculateThemeValue ( 4, tf[tfIndex][3], 2 );
		// Is the chart emphasis on water, cardinal, square aspects?
		if ( tfDominant[1] == 3 )		// water dominant?
			theme[3] += 1;
		if ( tfDominant[2] == 0 )		// cardinal dominant?
			theme[3] += 1;
		if ( dominantAspect == 3 )		// squares dominant
		theme[3] += 1;
		// ##### end Theme 4	#####
	
		// ##### Theme 5 #####
		calculateThemeValue ( 5, tf[tfIndex][4], 2 );
		// Is the chart emphasis on fire, fixed, trine aspects?
		if ( tfDominant[1] == 0 )		// fire dominant?
			theme[4] += 1;
		if ( tfDominant[2] == 1 )		// fixed dominant?
			theme[4] += 1;
		if ( dominantAspect == 2 )		// trines dominant
			theme[4] += 1;
		// ##### end Theme 5	#####
	
		// ##### Theme 6 #####
		calculateThemeValue ( 6, tf[tfIndex][5], 1 );
		// Is the chart emphasis on earth, mutable?
		if ( tfDominant[1] == 1 )		// earth dominant?
			theme[5] += 1;
		if ( tfDominant[2] == 2 )		// mutable dominant?
			theme[5] += 1;
		theme[5] = ( theme[5] > 3 ? theme[5]+1 : theme[5] );
		// ##### end Theme 6	#####
	
		// ##### Theme 7 #####
		calculateThemeValue ( 7, tf[tfIndex][6], 1 );
		// Is the chart emphasis on any of the following:  air, cardinal,
		// oppositions?
		if ( tfDominant[1] == 2 )		// air dominant?
			theme[6] += 1;
		if ( tfDominant[2] == 0 )		// cardinal dominant?
			theme[6] += 1;
		if ( dominantAspect == 1 )		// oppositions dominant
			theme[6] += 1;
		// ##### end Theme 7 #####
	
		// ##### Theme 8 #####
		var ruler = tf[tfIndex][7];
		if ( ruler == 9 )	// Pluto, modern ruler of Scorpio
		{
			calculateThemeValue ( 8, ruler, 0.5 );
			calculateThemeValue ( 8, 4, 0.5 );	// add contribution from ancient ruler Mars
			// also check Mars/Pluto aspects
			if ( orbType == 0 )
			{
				for ( n = 0; n < 7; n++ )
				{
					if ( isAspect ( planet[4], planet[ruler], a[n], ao[aoIndex][n] ))	// signRuler/Mars aspect
						theme[7] += ps[4]*ps[ruler]*aspectStrength ( planet[4], planet[ruler], a[n], ao[aoIndex][n], af[n] );
				}
			}
			else
			{
				orbValue = 0.5*(po[poIndex][4]+po[poIndex][ruler]);
				if ( isAspect ( planet[4], planet[ruler], a[n], orbValue ))	// signRuler/Mars aspect
					theme[7] += ps[4]*ps[ruler]*aspectStrength ( planet[4], planet[ruler], a[n], orbValue, af[n] );
			}
		}
		else
			calculateThemeValue ( 8, ruler, 1 );	// just use ancient ruler Mars
	
		// Is the chart emphasis on any of the following:  water, fixed
		if ( tfDominant[1] == 3 )		// water dominant?
			theme[7] += 1;
		if ( tfDominant[2] == 1 )		// fixed dominant?
			theme[7] += 1;
		theme[7] = ( theme[7] > 3 ? theme[7]+1 : theme[7] );
		// ##### end Theme 8 #####
	
		// ##### Theme 9 #####
		calculateThemeValue ( 9, tf[tfIndex][8], 1 );
		// Is the chart emphasis on any of the following:  fire, mutable, trine aspects
		if ( tfDominant[1] == 0 )		// fire dominant?
			theme[8] += 1;
		if ( tfDominant[2] == 2 )		// mutable dominant?
			theme[8] += 1;
		if ( dominantAspect == 2 )		// trines dominant
			theme[8] += 1;
		// ##### end Theme 9 #####
	
		// ##### Theme 10 #####
		calculateThemeValue ( 10, tf[tfIndex][9], 1 );
		// Is the chart emphasis on any of the following:  fire, mutable, trine aspects
		if ( tfDominant[1] == 1 )		//  earth dominant?
			theme[9] += 1;
		if ( tfDominant[2] == 0 )		// cardinal dominant?
			theme[9] += 1;
		if ( dominantAspect == 3 )		// squares dominant
			theme[9] += 1;
		// ##### end Theme 10 #####
	
		// ##### Theme 11 #####
		var ruler = tf[tfIndex][10];
		if ( ruler == 7 )	// Uranus, modern ruler of Aquarius
		{
			calculateThemeValue ( 11, ruler, 0.5 );
			calculateThemeValue ( 11, 6, 0.5 );	// add contribution from ancient ruler Saturn
			// also check Saturn/Uranus aspects
			if ( orbType == 0 )
			{
				for ( n = 0; n < 7; n++ )
				{
					if ( isAspect ( planet[6], planet[ruler], a[n], ao[aoIndex][n] ))	// signRuler/Mars aspect
						theme[10] += ps[6]*ps[ruler]*aspectStrength ( planet[6], planet[ruler], a[n], ao[aoIndex][n], af[n] );
				}
			}
			else
			{
				orbValue = 0.5*(po[poIndex][6]+po[poIndex][ruler]);
				if ( isAspect ( planet[6], planet[ruler], a[n], orbValue ))	// signRuler/Mars aspect
					theme[10] += ps[6]*ps[ruler]*aspectStrength ( planet[6], planet[ruler], a[n], orbValue, af[n] );
			}
		}
		else
			calculateThemeValue ( 11, ruler, 1 );	// just use ancient ruler Saturn
		// Is the chart emphasis on air, fixed, sextile aspects?
		if ( tfDominant[1] == 2 )		// air dominant?
			theme[10] += 1;
		if ( tfDominant[2] == 1 )		// fixed dominant?
			theme[10] += 1;
		if ( dominantAspect == 4 )		// sextiles dominant
			theme[10] += 1;
		// ##### end Theme 11 #####
	
		// ##### Theme 12 #####
		var ruler = tf[tfIndex][11];
		if ( ruler == 8 )	// Neptune, modern ruler of Pisces
		{
			calculateThemeValue ( 12, ruler, 0.5 );
			calculateThemeValue ( 12, 5, 0.5 );	// add contribution from ancient ruler Jupiter
			// also check Jupiter/Neptune aspects
			if ( orbType == 0 )
			{
				for ( n = 0; n < 7; n++ )
				{
					if ( isAspect ( planet[5], planet[ruler], a[n], ao[aoIndex][n] ))	// signRuler/Mars aspect
						theme[11] += ps[5]*ps[ruler]*aspectStrength ( planet[5], planet[ruler], a[n], ao[aoIndex][n], af[n] );
				}
			}
			else
			{
				orbValue = 0.5*(po[poIndex][5]+po[poIndex][ruler]);
				if ( isAspect ( planet[5], planet[ruler], a[n], orbValue ))	// signRuler/Mars aspect
					theme[11] += ps[5]*ps[ruler]*aspectStrength ( planet[5], planet[ruler], a[n], orbValue, af[n] );
			}

		}
		else
			calculateThemeValue ( 12, ruler, 1 );	// just use ancient ruler Saturn
		// Is the chart emphasis on water, mutable?
		if ( tfDominant[1] == 3 )		// water dominant?
			theme[11] += 1;
		if ( tfDominant[2] == 2 )		// mutable dominant?
			theme[11] += 1;
		theme[11] = ( theme[11] > 3 ? theme[11]+1 : theme[11] );
		// ##### end Theme 12 #####
	
		// add balance of polarities to each theme, 1 extra point distributed over
		// relevant themes
		if ( tfDominant[0] == 1 )
		{
			n = 0;
			while ( n < 12 )
			{
				theme[n] += 1/12;
				n += 2;
			}
		}
		else
		{
			n = 1;
			while ( n < 12 )
			{
				theme[n] += 1/12;
				n += 2;
			}
		}
		
		if ( precessionFlag  != 0 )
		{	// the First Point of Aries is not well-defined
			// adjusted var name to "nativityYear"
			m = Math.round(precession ( nativityYear  ) / 360  + 0.5 ); 	// also flaky
			for ( n = 0; n < 12; n++ )
			{
				k = ( n-m < 0 ? 12-m : n-m );
				precessedTheme[k] = theme[n];
			}
		}

		if ( chartType == 1 )	// moving 3-theme average: -1, 0, 1
		{	// display purposes only - DON'T HIT 'SAVE' 	!!!
			// find 3-theme moving average 4 times
			// quadruplicities * triplicities
			for ( m = 0; m < 4; m++ )
			{
				movingAverage ( theme );
				for ( n= 0 ; n < 12; n++ )
					theme[n] = avgThemeVal[n];
			}
		}
	}

	method = 1;	// specify default profile method

	function xProfile ( sArrayT, tArrayT )
	{
		var i; var  j; var  k;
		var n;
		var pCoincidence;
		// copies not pointers!
		// Type0: copy of raw data
		var sArrayType0 = [0,0,0,0,0,0,0,0,0,0,0,0];
		var tArrayType0 = [0,0,0,0,0,0,0,0,0,0,0,0];
		// Type1: normalised moving-average data
		var sArrayType1 = [[0,1,2,3,4,5,6,7,8,9,10,11],[0,0,0,0,0,0,0,0,0,0,0,0]];
		var tArrayType1 = [[0,1,2,3,4,5,6,7,8,9,10,11],[0,0,0,0,0,0,0,0,0,0,0,0]];
		// theme truth tables
		var sThemes = [0,0,0,0,0,0,0,0,0,0,0,0];
		var tThemes = [0,0,0,0,0,0,0,0,0,0,0,0];

		for ( n = 0; n < 12; n++ )
		{
			sArrayType0[n] = sArrayT[1][n];
			tArrayType0[n] = TidyUpAndFloat(tArrayT[1][n]);
			sArrayType1[1][n] = sArrayT[1][n];
			tArrayType1[1][n] = TidyUpAndFloat(tArrayT[1][n]);
		}
		// theme min and max values for (S)ubject and (T)arget from raw (type 0 ) data
		var tMinS;
		var tMinT;
		var tMaxS;
		var tMaxT;
		tMinS = locateThemeMin ( sArrayType0, 0, 12 );
		tMinT = locateThemeMin ( tArrayType0,0, 12 );
		tMaxS = locateThemeMax ( sArrayType0, 0, 12 );
		tMaxT = locateThemeMax ( tArrayType0, 0, 12 );
//debugP("S min "+tMinS+" max "+tMaxS);
//debugP("T min "+tMinT+" max "+tMaxT);
		// find closest correspondence between min and max in natal chart data (type 0)
		pCoincidence = minOf ( Math.abs ( tMinS - tMaxT ), Math.abs ( tMaxS - tMinT ) );
// use av of relative p/t separations? No, this queers the tPP/pCoincidence test below
//pCoincidence = ( Math.abs ( tMinS - tMaxT ) + Math.abs ( tMaxS - tMinT ) ) * 0.5;
		var tPP; var  tTT;	// diff between peaks
		tPP = Math.abs ( tMaxS - tMaxT );
		tTT = Math.abs ( tMinS - tMinT );	// currently unused
		// horror! suppose the PP coincidence is 0!
		if ( tPP < pCoincidence )
			pCoincidence = - ( tPP+1 );	// force lower range limit to -1 to distinguish from P/T 0

		// subject / target natal statistical data returned by statAnalyse
		var sStatsNatal = { mean: 0, sd: 0, skew: 0 };
		var tStatsNatal = { mean: 0, sd: 0, skew: 0 };
		statAnalyse ( sArrayType0, 0, 12, sStatsNatal );
		statAnalyse ( tArrayType0, 0, 12, tStatsNatal );
//debugP("natal skew: s = "+pNum(sStatsNatal.skew, precision)+" t = "+pNum(tStatsNatal.skew, precision));

		// find 3-theme moving average 4 times
		// quadruplicities * triplicities
		for ( m = 0; m < 4; m++ )
		{
			movingAverage ( sArrayType1[1] );
			for ( n= 0 ; n < 12; n++ )
				sArrayType1[1][n] = avgThemeVal[n];
		
			movingAverage ( tArrayType1[1] );
			for ( n= 0 ; n < 12; n++ )
				tArrayType1[1][n] = avgThemeVal[n];
		}
		// normalise
		var themeMax = 0;
		for ( n = 0; n < 12; n++ )
			if ( sArrayType1[1][n] > themeMax )
				themeMax = sArrayType1[1][n];

		for ( n = 0; n < 12; n++ )
			sArrayType1[1][n] /= themeMax;

		themeMax = 0;
		for ( n = 0; n < 12; n++ )
			if ( tArrayType1[1][n] > themeMax )
			themeMax = tArrayType1[1][n];

		for ( n = 0; n < 12; n++ )
			tArrayType1[1][n] /= themeMax;

		var data = { fitC: 0, fitS: 0, pCoincidence: 0 };	// pre-load meaningless values
		data.pCoincidence = pCoincidence;	// natal correspondence
		// we DO NOT deal with trough correspondence since the peak values are by definition primary
		// if the two peaks _or_ troughs are in close correspondence this discourages close relationship
		// well, this seems so for P-P correspondence
/*
{
debugP("mean: s = "+pNum(sStats.mean, precision)+" t = "+pNum(tStats.mean, precision));
debugP("sigma: s = "+pNum(sStats.sd, precision)+" t = "+pNum(tStats.sd, precision));
debugP("mean-/+sigma: s = "+pNum((sStats.mean-sStats.sd), precision)+", "+pNum((sStats.mean+sStats.sd), precision)+" t = "+pNum((tStats.mean-tStats.sd), precision)+", "+pNum((tStats.mean+tStats.sd), precision));
debugP("rolloff rate: s = "+pNum((sStats.sd/sStats.mean), precision)+" t = "+pNum((tStats.sd/tStats.mean), precision));
		// rate at which sigma declines from mean
		// larger rolloff means tighter peak
		// larger skew means greater/fatter tail
		// skew negative means L tail dominant, else R tail
debugP("skew: s = "+pNum(sStats.skew, precision)+" t = "+pNum(tStats.skew, precision));
}
*/
//debugP("transformed skew: s = "+pNum(sStats.skew, precision)+" t = "+pNum(tStats.skew, precision));

		if ( method == 1 )
		{
			// inflection point value, location
			var sArrayTp = [[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0]];
			var tArrayTp = [[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0]];
			// (C)omplementarity / (S)imilarity of subject / target Type1 data
			var dataL;		// linear,eq expression coefficients
			var sSlope; var  tSlope;	// subject/target dataL.slope
			// initialise num turning points found
			var numInflectionsS = 0;
			var numInflectionsT = 0;
			var lastSlopeS; var  lastSlopeT;
			// stat package here is for broad relationship analysis ignoring minor S/T differences
			var sStats = { mean: 0, sd: 0, skew: 0 };
			var tStats = { mean: 0, sd: 0, skew: 0 };
			// standard statistical analysis of transformed data
			statAnalyse ( sArrayType1[1], 0, 12, sStats );
			statAnalyse ( tArrayType1[1], 0, 12, tStats );
			// data used by themeFit determining bounds of insignificance
			sStats.sd = sStats.sd*sdFactor;
			tStats.sd = tStats.sd*sdFactor;

			tMinS1 = locateThemeMin ( sArrayType1[1], 0, 12 );
			tMinT1 = locateThemeMin ( tArrayType1[1],0, 12 );
			tMaxS1 = locateThemeMax ( sArrayType1[1], 0, 12 );
			tMaxT1 = locateThemeMax ( tArrayType1[1], 0, 12 );
		
			// initialise slopes to ends of graph
			if ( sArrayType1[1][11] > sArrayType1[1][0] )
				lastSlopeS = -1;
			else
				if ( sArrayType1[1][11] < sArrayType1[1][0] )
					lastSlopeS = 1;
				else
					lastSlopeS = 0;
		
			if ( tArrayType1[1][11] > tArrayType1[1][0] )
				lastSlopeT = -1;
			else
				if ( tArrayType1[1][11] < tArrayType1[1][0] )
					lastSlopeT = 1;
				else
					lastSlopeT = 0;
		
			for ( i = 0; i < 11; i++ )
			{
				// with only 2 values this is excessive
				dataL = fitLinearEq ( sArrayType1[0], sArrayType1[1], i, 2 );
				sSlope = dataL.slope;
				dataL = fitLinearEq ( tArrayType1[0], tArrayType1[1], i, 2 );
				tSlope = dataL.slope;

				if ( lastSlopeS != signOf(sSlope) )
				{
					sArrayTp[0][numInflectionsS]= i;
					sArrayTp[1][numInflectionsS] = signOf(sSlope);
					numInflectionsS++;
				}
		
				if ( lastSlopeT != signOf(tSlope) )
				{
					tArrayTp[0][numInflectionsT] = i;
					tArrayTp[1][numInflectionsT]= signOf(tSlope);
					numInflectionsT++;
				}
				lastSlopeS = signOf(sSlope);
				lastSlopeT = signOf(tSlope);
			}

			// findLinear doesn't	 wrap so check array end against start here
			if ( sArrayType1[1][11] > sArrayType1[1][0] )
				sSlope = -1;
			else
				if ( sArrayType1[1][11] < sArrayType1[1][0] )
					sSlope = 1;
				else
					sSlope = 0;
		
			if ( tArrayType1[1][11] > tArrayType1[1][0] )
				tSlope = -1;
			else
				if ( tArrayType1[1][11] < tArrayType1[1][0] )
					tSlope = 1;
				else
					tSlope = 0;

			if ( ( lastSlopeS != sSlope ) && ( sSlope != 0 ) )
			{
				sArrayTp[0][numInflectionsS]= 11;
				sArrayTp[1][numInflectionsS] = sSlope;
			}
	
			if ( ( lastSlopeT != tSlope ) && ( tSlope != 0 ) )
			{
				tArrayTp[0][numInflectionsT] = 11;
				tArrayTp[1][numInflectionsT]= tSlope;
			}
// these need exporting to returned method 2 package!
// this will only happen if we call method 1 then method 2 since the data isn't preserved
// between s/t pair calls
			var scale;

			for ( i = 0; i < numInflectionsS; i++ )
			{
				// if turning point +ve find trough limits else find peak limits; ignore slope 0 regions
				if ( sArrayTp[1][i] == 1 )	// find trough
				{
					troughLimits ( sArrayType1[1], sArrayTp[0][i], sThemes );

					for ( j = 0; j < numInflectionsT; j++ )
					{
						scale = scaleFactor ( sArrayTp[0][i], tArrayTp[0][j] );

						if ( tArrayTp[1][j] == 1 )	// find trough
						{
							troughLimits ( tArrayType1[1], tArrayTp[0][j], tThemes );

							for ( k = 0; k < 12; k++ )
								if ( sThemes[k] & tThemes[k] )
									data.fitS += scale * themeFit ( sArrayType1[1][k], tArrayType1[1][k], sArrayTp[1][i], tArrayTp[1][j], sStats, tStats );
						}
						else
						{
							if ( tArrayTp[1][j] == -1 )	// find peak
							{
								peakLimits ( tArrayType1[1], tArrayTp[0][j], tThemes );

								for ( k = 0; k < 12; k++ )
									if ( sThemes[k] & tThemes[k] )
										data.fitC += scale * themeFit ( sArrayType1[1][k], tArrayType1[1][k], sArrayTp[1][i], tArrayTp[1][j], sStats, tStats );
							}
						}
					}
				}
				else
				{
					if ( sArrayTp[1][i] == -1 )	// find peak
					{
						peakLimits ( sArrayType1[1], sArrayTp[0][i], sThemes );
			
						for ( j = 0; j < numInflectionsT; j++ )
						{
							scale = scaleFactor ( sArrayTp[0][i], tArrayTp[0][j] );

							if ( tArrayTp[1][j] == 1 )	// find trough
							{
								troughLimits ( tArrayType1[1], tArrayTp[0][j], tThemes );

								for ( k = 0; k < 12; k++ )
									if ( sThemes[k] & tThemes[k] )
										data.fitC += scale * themeFit ( sArrayType1[1][k], tArrayType1[1][k], sArrayTp[1][i], tArrayTp[1][j], sStats, tStats	);
							}
							else
							{
								if ( tArrayTp[1][j] == -1 )	// find peak
								{
									peakLimits ( tArrayType1[1], tArrayTp[0][j], tThemes );

									for ( k = 0; k < 12; k++ )
										if ( sThemes[k] & tThemes[k] )
											data.fitS += scale * themeFit ( sArrayType1[1][k], tArrayType1[1][k], sArrayTp[1][i], tArrayTp[1][j], sStats, tStats	);
								}
							}
						}
					}
				}
			}

			if ( ( data.fitC == 0 ) && ( data.fitS == 0 ) )	// fails if both are exactly zero
				debugP("Error: sdFactor too large, no significant data found");

debugP("fitC "+pNum(data.fitC, precision)+" fitS "+pNum(data.fitS, precision)+" pCoincidence "+data.pCoincidence);
			return data;	// struct {fitC, fitS, pCoincidence }
		}	// end method 1
		
		if ( method == 2 )
		{
			// dummy data till I get something sensible returned curveData et.al.)
			var data = { fitC: 0, fitS: 0, pCoincidence: 0 };	// dummy for fn return
			// test run on Type1 data - find primary peaks
			// standard statistical analysis of transformed data
			var sStats = { mean: 0, sd: 0, skew: 0 };
			var tStats = { mean: 0, sd: 0, skew: 0 };
			// standard statistical analysis of transformed data
			statAnalyse ( sArrayType1[1], 0, 12, sStats );
			statAnalyse ( tArrayType1[1], 0, 12, tStats );
			// data used by themeFit determining bounds of insignificance
			sStats.sd = sStats.sd*sdFactor;
			tStats.sd = tStats.sd*sdFactor;
			// find closest correspondence between min and max
			tMinS2 = locateThemeMin ( sArrayType1[1], 0, 12 );
			tMinT2 = locateThemeMin ( tArrayType1[1],0, 12 );
			tMaxS2 = locateThemeMax ( sArrayType1[1], 0, 12 );
			tMaxT2 = locateThemeMax ( tArrayType1[1], 0, 12 );

var skewRatio = sStats.skew / tStats.skew;

			// note the derived coefficients only apply over the range of values analysed
			// centre values are min/max values (above)
//debugP("Subject peak");
			limits = peakLimits ( sArrayType1[1], tMaxS2, sThemes );
//debugP(sThemes);
			sDataP = curveAnalyse ( sArrayType1[0], sArrayType1[1], tMaxS2, limits );
//predictValues ( limits, sDataP, sArrayType1[1] )

//debugP("Subject trough");
			limits = troughLimits ( sArrayType1[1], tMinS2, sThemes );
//debugP(sThemes);
			sDataT = curveAnalyse ( sArrayType1[0], sArrayType1[1], tMinS2, limits );
//predictValues ( limits, sDataT, sArrayType1[1] )

//debugP("Target peak");
			limits = peakLimits ( tArrayType1[1], tMaxT2, tThemes );
//debugP(tThemes);
			tDataP = curveAnalyse ( tArrayType1[0], tArrayType1[1], tMaxT2, limits );
//predictValues ( limits, tDataP, tArrayType1[1] )

//debugP("Target trough");
			limits = troughLimits ( tArrayType1[1], tMinT2, tThemes );
//debugP(tThemes);
			tDataT = curveAnalyse ( tArrayType1[0], tArrayType1[1], tMinT2, limits );
//predictValues ( limits, tDataT, tArrayType1[1] )
debugP("skew "+pNum(sStats.skew, precision)+" "+pNum(tStats.skew, precision));
/*
			var rejected = 0;
			//		if (  minOf ( Math.abs ( sDataP.centre - tDataT.centre ), Math.abs ( sDataT.centre - tDataP.centre ) )  > 1 )
			rejected =  (  minOf ( Math.abs ( sDataP.centre - tDataT.centre ), Math.abs ( sDataT.centre - tDataP.centre ) )  > 1 ? 1 : 0 );

			if ( rejected )
debugP("Rejected: lacks close correspondence in P/T" );
//if ( minOf ( Math.abs ( sDataP.centre - tDataP.centre ),  Math.abs ( sDataT.centre - tDataT.centre ) < 4 ) )
			if (!rejected)
			rejected =  ( Math.abs ( sDataP.centre - tDataP.centre ) < 2 ? 1 : 0 );

			if ( rejected )
			{
// debugP("Rejected, P/P or T/T correspondence, difficult");
debugP("Rejected, P/P correspondence, difficult");
//debugP("min PP "+Math.abs ( sDataP.centre - tDataP.centre ));
// debugP("min TT "+Math.abs ( sDataT.centre - tDataT.centre ));
			}
// even given these are satisfied, have to check relative weighting of of non-P/T areas.
// if the difference is large, this mitigates against closenes as implies complementarity
			if (!rejected )
			rejected = ( signOf ( sStats.skew ) != signOf ( tStats.skew ) ? 1 : 0 );
						
			if ( rejected )
				debugP("Rejected: weighting complementary, friend not partner");
// fails on rejected - only gets last test result...

			if ( !rejected ) // use stats approach				
			{
//debugP("skewRatio "+pNum(skewRatio, precision));
debugP("Possible...");
//					if ( signOf ( sStats.skew ) != signOf ( tStats.skew ) )
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
//				}
			}
*/
// minOf here just gives the minimum theme not min distance!
// should we be referring to pCoincidence AND p centres?
//			if ( ( pCoincidence < 3 ) & ( signOf ( sStats.skew ) == signOf ( tStats.skew ) ) & ( Math.abs ( sDataP.centre - tDataP.centre ) > 1 ) )
// note we check skew of averaged data here since we are interested in major features only

if ( ( pCoincidence < 3 ) && ( signOf ( sStats.skew ) == signOf ( tStats.skew ) ) )
			{
if ( signOf ( sStats.skew ) != signOf ( sStatsNatal.skew ) )
	debugP ("Warning: subject natal skew changed from "+pNum(sStatsNatal.skew, precision)+" to "+pNum(sStats.skew, precision) );
if ( signOf ( tStats.skew ) != signOf ( tStatsNatal.skew ) )
	debugP ("Warning: target natal skew changed from "+pNum(tStatsNatal.skew, precision)+" to "+pNum(tStats.skew, precision) );
				debugP("passes natal data check");
				if ( ( pCoincidence < 0 ) && ( pCoincidence > -4 ) )
					debugP("but natal peaks too close, separation "+(Math.abs(pCoincidence+1))+", failed");
			}			
			else
				debugP("Failed");

			return data;	// send default rubbish back for record.js or profile popup auto-closes
		}	// end method 2
	}

// END OF ENGINE
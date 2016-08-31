var NITFHeader = "";

function readBlob(opt_startByte, opt_stopByte) {
	var files = document.getElementById('files').files;
	var file = files[0];
	var start = parseInt(opt_startByte) || 0;
	var stop = parseInt(opt_stopByte) || file.size - 1;
	var reader = new FileReader();
	var text;
	// If we use onloadend, we need to check the readyState.
	reader.onloadend = function(evt) {
		if (evt.target.readyState == FileReader.DONE) { // DONE == 2

			NITFHeader = evt.target.result;
			checkFormat();
		}
	};
	var blob = file.slice(start, stop + 1);
	reader.readAsText(blob);


}



function readNITFMetadata(){
	readBlob(0, 3600); //read first 3600 bytes in file
}

function checkFormat(){

	if (NITFHeader.indexOf("CMETAA") >= 0) {
		parseNITF(NITFHeader);
	}else if (NITFHeader.indexOf("SICD") >= 0) {
		parseSICD(NITFHeader);
	}else{
		alert('Selected file is not a complex NITF or SICD');
	}
	
}

function parseSICD(NITFHeader){

	var headerLength = parseInt(NITFHeader.substr(354, 6)); // File header length
	var numI = parseInt(NITFHeader.substr(360, 3)); // Number of image segments

	var byteCounter = 363; //need to keep track of byte pointer to make sure we find the right byte of the different offsets
	var ImgSubhdrLngths = 0;
	var ImgLngths = 0;
	if (numI > 0){
		for (var i=0;i<numI;i++){
			ImgSubhdrLngths = ImgSubhdrLngths + parseInt(NITFHeader.substr(byteCounter, 6));
			byteCounter = byteCounter + 6;
			ImgLngths = ImgLngths + parseInt(NITFHeader.substr(byteCounter, 10));
			byteCounter = byteCounter + 10;
		}
	}

	// Jump over segments not used by SICD to DES description
	byteCounter = byteCounter + 6;

	// Jump over text segments
	var NUMT = parseInt(NITFHeader.substr(byteCounter, 3));
	byteCounter = byteCounter + 3;
	
	var LTSH = 0;
	var LT = 0;
	for (var lp=0; lp<NUMT; lp++){
		LTSH = LTSH + parseInt(NITFHeader.substr(byteCounter, 4)); // Length of text segment subheader
		byteCounter = byteCounter + 4;
		LT = LT + parseInt(NITFHeader.substr(byteCounter, 5)); // Length of text segment data
		byteCounter = byteCounter + 5;
	}
	var NUMDES = parseInt(NITFHeader.substr(byteCounter, 3));  // Number of data extension segments
	byteCounter = byteCounter + 3;
	
	LDSH = [];
	LD = [];
	for (lp=0; lp<NUMDES; lp++){
		LDSH[lp] = parseInt(NITFHeader.substr(byteCounter, 4)); // Length of data extension segment subheader
		byteCounter = byteCounter + 4;
		LD[lp] = parseInt(NITFHeader.substr(byteCounter, 9)); // Length of data extension segment data
		byteCounter = byteCounter + 9;
	}
	// SICD Volume 2, File Format Description, section 3.1.1 says that SICD XML
	// metadata must be stored in first DES.
	var SICD_DES_Offset = headerLength + ImgSubhdrLngths + ImgLngths + LTSH + LT + LDSH[0];
	var SICD_DES_Length = LD[0];


	readSICD_XML(SICD_DES_Offset, SICD_DES_Offset + SICD_DES_Length);
	
	
}



function readSICD_XML(opt_startByte, opt_stopByte) {
	var files = document.getElementById('files').files;
	var file = files[0];
	var start = parseInt(opt_startByte) || 0;
	var stop = parseInt(opt_stopByte) || file.size - 1;
	var reader = new FileReader();
	var text;
	// If we use onloadend, we need to check the readyState.
	reader.onloadend = function(evt) {
		if (evt.target.readyState == FileReader.DONE) { // DONE == 2
			parseSICD_XML(evt.target.result);
		}
	};
	var blob = file.slice(start, stop + 1);
	reader.readAsText(blob);

}


function parseSICD_XML(SICD_XML){
	
	// look for SICD string in XML
	if (SICD_XML.substr(0,5) == "<SICD") {
		SICD_XML = SICD_XML.replace(/>\s*/g,">");
		SICD_XML = SICD_XML.replace(/\s*</g,"<");
		//parse out the values from the XML
		parser= new DOMParser();
		xmlDoc = parser.parseFromString(SICD_XML, "text/xml");
		
		
		var IID = xmlDoc.getElementsByTagName("CoreName")[0].childNodes[0].nodeValue;

		var CollectStart = xmlDoc.getElementsByTagName("CollectStart")[0].childNodes[0].nodeValue;
		
		var year = parseInt(CollectStart.substr(0,4));
		var month = parseInt(CollectStart.substr(5,2));
		var day = parseInt(CollectStart.substr(8,2));
		var hour = parseInt(CollectStart.substr(11,2));
		var minute = parseInt(CollectStart.substr(14,2));
		var second = parseInt(CollectStart.substr(17,2));
		var ms = parseInt(CollectStart.substr(20,3));
		
		var imageDateTime = Date.UTC(year, month - 1, day, hour, minute, second, ms); //milliseconds since Jan 1, 1970 in UTC (month starts at 0, not 1)
		
		var temp = new Date(imageDateTime).toUTCString();
				
		var collectDuration = parseFloat(xmlDoc.getElementsByTagName("CollectDuration")[0].childNodes[0].nodeValue);
		var ARPPolyX = [];
		var ARPPolyY = [];
		var ARPPolyZ = [];

		var numXCoeffs = parseInt(xmlDoc.getElementsByTagName("ARPPoly")[0].childNodes[0].attributes[0].nodeValue);
		var numYCoeffs = parseInt(xmlDoc.getElementsByTagName("ARPPoly")[0].childNodes[1].attributes[0].nodeValue);
		var numZCoeffs = parseInt(xmlDoc.getElementsByTagName("ARPPoly")[0].childNodes[2].attributes[0].nodeValue);
		for( var i=0; i<numXCoeffs+1; i++){
			ARPPolyX[i] = parseFloat(xmlDoc.getElementsByTagName("ARPPoly")[0].childNodes[0].childNodes[i].textContent);
		}
		for( var i=0; i<numYCoeffs+1; i++){
			ARPPolyY[i] = parseFloat(xmlDoc.getElementsByTagName("ARPPoly")[0].childNodes[1].childNodes[i].textContent);
		}
		for( var i=0; i<numZCoeffs+1; i++){
			ARPPolyZ[i] = parseFloat(xmlDoc.getElementsByTagName("ARPPoly")[0].childNodes[2].childNodes[i].textContent);
		}
		var LatDD = parseFloat(xmlDoc.getElementsByTagName("LLH")[0].childNodes[0].childNodes[0].textContent);
		var LonDD = parseFloat(xmlDoc.getElementsByTagName("LLH")[0].childNodes[1].childNodes[0].textContent);
		var HAE = parseFloat(xmlDoc.getElementsByTagName("LLH")[0].childNodes[2].childNodes[0].textContent);
				
	
		processMetadata(HAE, LatDD, LonDD, collectDuration, ARPPolyX, ARPPolyY, ARPPolyZ, imageDateTime, IID);
		
		
	}else{
		alert('Image is not a SICD');
	}

}

function parseNITF(data){

	var rad2deg = 180/Math.PI;
	var deg2rad = Math.PI/180;

	var headerLength = data.substr(354, 6);
	
	var GHBlk30 = false;
	var GHBlk40 = false;
	var ASARS = false;
	var LynxSAR = false;
	
	
	var ACFTB = data.substr(data.indexOf("ACFTB")+11, parseInt(data.substr(data.indexOf("ACFTB")+6, 5)));
	var sensorType = ACFTB.substr(46,6).trim();

	if(sensorType == "RTNESS"){
		GHBlk30 = true;
	}else if(sensorType == "MPRTIP"){
		GHBlk40 = true;
	}else if(sensorType == "AIP"){
		ASARS = true;
	}else if(sensorType == "ANAPY8"){
		LynxSAR = true;
	}
	
	var CMETAA = data.substr(data.indexOf("CMETAA")+11, parseInt(data.substr(data.indexOf("CMETAA")+6, 5)));	
	var TRETag = "IMComposite";
	if (data.indexOf("IMComposite") == -1) {
		TRETag = "IMNot Avail"; //LynxSAR TRE Tag is defined as IMNotAvail instead of IMComposite
	}
	
	var IID = data.substr(data.indexOf(TRETag)+43, 16);
	
	var collectionTime = data.substr(data.indexOf(TRETag)+12, 14);
	var year = parseInt(collectionTime.substr(0,4));
	var month = parseInt(collectionTime.substr(4,2));
	var day = parseInt(collectionTime.substr(6,2));
	var hour = parseInt(collectionTime.substr(8,2));
	var minute = parseInt(collectionTime.substr(10,2));
	var second = parseInt(collectionTime.substr(12,2));
	
	var imageDateTime = Date.UTC(year, month - 1, day, hour, minute, second); //milliseconds since Jan 1, 1970 in UTC (month starts at 0, not 1)
	
	var coordType = CMETAA.substr(847,5).trim();
	var collectDuration = parseFloat(CMETAA.substr(1555,7));
	
	
	var APCEN_X = parseFloat(CMETAA.substr(903,13));
	var APCEN_Y = parseFloat(CMETAA.substr(916,13));
	var APCEN_Z = parseFloat(CMETAA.substr(929,13));
	
	if (coordType == "WGS84"){
		var LatDD = parseFloat(CMETAA.substr(1008,13));
		var LonDD = parseFloat(CMETAA.substr(1021,13));
		var HAE = parseFloat(CMETAA.substr(1034,13));

		//convert everything from WGS84 to ECEF
		var e2 = .0066943799901377997; // eccentricity squared of Earth (WGS 84 value)
		var a = 6378137.0;              // semimajor radius of the Earth (WGS 84 value)

		// calculate distance to surface of ellipsoid
		var R = a / Math.sqrt(1.0 - e2 * Math.sin(APCEN_X*deg2rad) * Math.sin(APCEN_X*deg2rad));

		// calculate coordinates
		var ARPPos_X = (R + APCEN_Z) * Math.cos(APCEN_X*deg2rad) * Math.cos(APCEN_Y*deg2rad);
		var ARPPos_Y = (R + APCEN_Z) * Math.cos(APCEN_X*deg2rad) * Math.sin(APCEN_Y*deg2rad);
		var ARPPos_Z = (R + APCEN_Z - e2 * R) * Math.sin(APCEN_X*deg2rad);
		
	}else if(coordType == "ECEF"){
		
		var ECEF_X = parseFloat(CMETAA.substr(1008,13));
		var ECEF_Y  = parseFloat(CMETAA.substr(1021,13));
		var ECEF_Z  = parseFloat(CMETAA.substr(1034,13));
		var LLA = ECEFtoLLA(ECEF_X, ECEF_Y, ECEF_Z);
		var LatDD = LLA.lat;
		var LonDD = LLA.lon;
		var HAE = LLA.alt;

		var ARPPos_X = APCEN_X;
		var ARPPos_Y = APCEN_Y;
		var ARPPos_Z = APCEN_Z;
	}

	

	var ARPVel_X = parseFloat(CMETAA.substr(1067,10));
	var ARPVel_Y = parseFloat(CMETAA.substr(1077,10));
	var ARPVel_Z = parseFloat(CMETAA.substr(1087,10));
	var SNACC_X = parseFloat(CMETAA.substr(1097,10));
	var SNACC_Y = parseFloat(CMETAA.substr(1107,10));
	var SNACC_Z = parseFloat(CMETAA.substr(1117,10));

	/// ECEF/NED adjustments
	// Some fields in the NITF TREs (like velocity) appear to be
	// north-east-down, not ECEF so convert NED to ECEF
	if (ASARS){
		//scene center
		var SCECN_X = parseFloat(CMETAA.substr(1008,13));
		var SCECN_Y = parseFloat(CMETAA.substr(1021,13));
		var SCECN_Z = parseFloat(CMETAA.substr(1034,13));
		
		//convert NED to ECEF
		var orp_ecf = [ SCECN_X, SCECN_Y, SCECN_Z ];
		var vel_ned = [ ARPVel_X, ARPVel_Y, ARPVel_Z ];
		var acc_ned = [ SNACC_X, SNACC_Y, SNACC_Z ];
		var vel_ecf = ned_to_ecf(vel_ned, orp_ecf);
		var acc_ecf = ned_to_ecf(acc_ned,orp_ecf);
		
		//save adjusted values
		ARPVel_X = vel_ecf[0];
		ARPVel_Y = vel_ecf[1];
		ARPVel_Z = vel_ecf[2];
		SNACC_X = acc_ecf[0];
		SNACC_Y = acc_ecf[1];
		SNACC_Z = acc_ecf[2];
		
	}

	// Move current position to start position since ARPPos is based on center aperture (AP_CEN)
	var SCPTime = collectDuration/2.0;
	var ARPPolyX = [ ARPPos_X - (ARPVel_X * SCPTime) + ((SNACC_X/2) * Math.pow(SCPTime,2)), ARPVel_X - (SNACC_X * SCPTime), SNACC_X/2, 0, 0, 0 ];
	var ARPPolyY = [ ARPPos_Y - (ARPVel_Y * SCPTime) + ((SNACC_Y/2) * Math.pow(SCPTime,2)), ARPVel_Y - (SNACC_Y * SCPTime), SNACC_Y/2, 0, 0, 0 ];
	var ARPPolyZ = [ ARPPos_Z - (ARPVel_Z * SCPTime) + ((SNACC_Z/2) * Math.pow(SCPTime,2)), ARPVel_Z - (SNACC_Z * SCPTime), SNACC_Z/2, 0, 0, 0 ];
	
	processMetadata(HAE, LatDD, LonDD, collectDuration, ARPPolyX, ARPPolyY, ARPPolyZ, imageDateTime, IID);

}

function processMetadata(HAE, LatDD, LonDD, CollectDuration, X, Y, Z, imageDateTime, IID){
	
	var earthEccentricitysqrd = 0.0066943799901377997;	//WGS-84 value
	var semiMajorRadius = 6378137.0;					//WGS-84 value
	var semiMinorRadius = 6356752.31424518;
	
	//convert to radians for calculations
	var Lat = LatDD * Math.PI/180;
	var Lon = LonDD * Math.PI/180;


	var R = semiMajorRadius/(Math.sqrt(1-(earthEccentricitysqrd*Math.sin(Lat)*Math.sin(Lat))));
	var test = (Math.sqrt(1-(earthEccentricitysqrd*Math.sin(Lat)*Math.sin(Lat))));
	var poi_x = (R+HAE)*Math.cos(Lat)*Math.cos(Lon);
	var poi_y = (R+HAE)*Math.cos(Lat)*Math.sin(Lon);
	var poi_z = (R+earthEccentricitysqrd-(R*earthEccentricitysqrd))*Math.sin(Lat);

	var normal_vector_x=poi_x/(semiMajorRadius*semiMajorRadius);
	var normal_vector_y=poi_y/(semiMajorRadius*semiMajorRadius);
	var normal_vector_z=poi_z/(semiMinorRadius*semiMinorRadius);

	var temp_magnitude = Math.sqrt(Math.pow(normal_vector_x, 2)+Math.pow(normal_vector_y, 2)+Math.pow(normal_vector_z, 2));

	var KDPx = normal_vector_x/temp_magnitude;
	var KDPy = normal_vector_y/temp_magnitude;
	var KDPz = normal_vector_z/temp_magnitude;

	temp_magnitude = Math.sqrt((Math.pow(KDPx, 2)+Math.pow(KDPy, 2)+Math.pow(KDPz, 2)));
	var nx = KDPx/temp_magnitude;
	var ny = KDPy/temp_magnitude;
	var nz = KDPz/temp_magnitude;

	var ProjNx = 0-(nz*nx);
	var ProjNy = 0-(nz*ny);
	var ProjNz = 1-(nz*nz);

	temp_magnitude = Math.sqrt((Math.pow(ProjNx, 2)+Math.pow(ProjNy, 2)+Math.pow(ProjNz, 2)));
	ProjNx = ProjNx/temp_magnitude;
	ProjNy = ProjNy/temp_magnitude;
	ProjNz = ProjNz/temp_magnitude;


	var azimuth_arr = new Array();
	var graze_arr = new Array();
	var i = 0;
	

	for (var time=0; time<CollectDuration; time=time+(CollectDuration/100)){

		var Px = X[5]*Math.pow(time, 5) + X[4]*Math.pow(time, 4) + X[3]*Math.pow(time, 3) + X[2]*Math.pow(time, 2) + X[1]*time + X[0];
		var Py = Y[5]*Math.pow(time, 5) + Y[4]*Math.pow(time, 4) + Y[3]*Math.pow(time, 3) + Y[2]*Math.pow(time, 2) + Y[1]*time + Y[0];
		var Pz = Z[5]*Math.pow(time, 5) + Z[4]*Math.pow(time, 4) + Z[3]*Math.pow(time, 3) + Z[2]*Math.pow(time, 2) + Z[1]*time + Z[0];
		
		var Rx = Px-poi_x;
		var Ry = Py-poi_y;
		var Rz = Pz-poi_z;

		temp_magnitude = Math.sqrt((Math.pow(Rx, 2)+Math.pow(Ry, 2)+Math.pow(Rz, 2)));
		Rx = Rx/temp_magnitude;
		Ry = Ry/temp_magnitude;
		Rz = Rz/temp_magnitude;

		var dotx = Rx - (((Rx*nx)+(Ry*ny)+(Rz*nz))*nx);
		var doty = Ry - (((Rx*nx)+(Ry*ny)+(Rz*nz))*ny);
		var dotz = Rz - (((Rx*nx)+(Ry*ny)+(Rz*nz))*nz);

		temp_magnitude = Math.sqrt((Math.pow(dotx, 2)+Math.pow(doty, 2)+Math.pow(dotz, 2)));
		var JDPx = dotx/temp_magnitude;
		var JDPy = doty/temp_magnitude;
		var JDPz = dotz/temp_magnitude;

		var dot_projN_jdp = (JDPx*ProjNx)+(JDPy*ProjNy)+(JDPz*ProjNz);
		var cross_jdp_projNx = (JDPy*ProjNz)-(JDPz*ProjNy);
		var cross_jdp_projNy = (JDPz*ProjNx)-(JDPx*ProjNz);
		var cross_jdp_projNz = (JDPx*ProjNy)-(JDPy*ProjNx);

		var dot_jdp_kdp = (cross_jdp_projNx*KDPx)+(cross_jdp_projNy*KDPy)+(cross_jdp_projNz*KDPz);

		var temp_az = 0;
		var temp_grz = 0;


		if (Math.atan2(dot_jdp_kdp, dot_projN_jdp)<0){
			temp_az = Math.atan2(dot_jdp_kdp,dot_projN_jdp)+(2*Math.PI);
		}else{
			temp_az = Math.atan2(dot_jdp_kdp,dot_projN_jdp);
		}

		temp_grz = Math.asin((Rx*KDPx)+(Ry*KDPy)+(Rz*KDPz));


		azimuth_arr[i] = Math.round(100.0*(temp_az * 180/Math.PI))/100.0;
		graze_arr[i] = Math.round(100.0*(temp_grz * 180/Math.PI))/100.0;

		i++;

	}


	var centerAzimuth = azimuth_arr[50];
	var centerGraze = graze_arr[50];

	var path_start_az = azimuth_arr[0];
	var path_middle_az = azimuth_arr[49];
	var path_end_az = azimuth_arr[99];
	
	var path_start_grz = graze_arr[0];
	var path_middle_grz = graze_arr[49];
	var path_end_grz = graze_arr[99];

	
	var imageData = {
			azimuth: centerAzimuth,
			collectDuration: CollectDuration,
			graze: centerGraze,
			collectionDate: imageDateTime,
			latDD: LatDD,
			lonDD: LonDD,
			path_start_az: path_start_az,
			path_middle_az: path_middle_az,
			path_end_az: path_end_az,
			path_start_grz: path_start_grz,
			path_middle_grz: path_middle_grz,
			path_end_grz: path_end_grz,
			imageId: IID
		};

	
}

function ECEFtoLLA(x, y, z){

		// define constants
		var e2 = 6.6943799901377997e-3;  // eccentricity squared of Earth (WGS 84 value)
		var a = 6378137.0;               // semimajor radius of the Earth (WGS 84 value)

		// calculate derived constants
		var e4 = e2 * e2;
		var ome2 = 1.0 - e2;
		var a2 = a * a;
		var b = a * Math.sqrt(ome2);         // semiminor radius of the Earth
		var b2 = b * b;
		var e_b2 = (a2 - b2) / b2;


		// calculate intermediates
		var z2 = z * z;
		var r2 = (x * x) + (y * y);
		var r = Math.sqrt(r2);

		// Check for invalid solution
		var valid = ((a * r) * (a * r) + (b * z) * (b * z) > (a2 - b2) * (a2 - b2));


		// calculate longitude
		var lon = Math.atan2(y, x)*180/Math.PI; // atan2d not available until MATLAB 2012b

		// calculate intermediates
		var F = 54.0 * b2 * z2;
		var G = r2 + ome2 * z2 - e2 * (a2 - b2);
		var c = e4 * F * r2 / (G * G * G);
		var s = Math.pow(1.0 + c + Math.sqrt(c * c + 2 * c), (1 / 3));
		var templ = s + 1.0 / s + 1.0;
		var P = F / (3.0 * templ * templ * G * G);
		var Q = Math.sqrt(1.0 + 2.0 * e4 * P);
		var r0 = -P * e2 * r / (1.0 + Q) + Math.sqrt(Math.abs(0.5 * a2 * (1.0 + 1.0 / Q) - P * ome2 * z2 / (Q * (1.0 + Q)) - 0.5 * P * r2));
		var temp2 = r - e2 * r0;
		var U = Math.sqrt(temp2 * temp2 + z2);
		var V = Math.sqrt(temp2 * temp2 + ome2 * z2);
		var z0 = b2 * z / (a * V);

		// calculate latitude
		var lat = Math.atan2(z + e_b2 * z0, r)*180/Math.PI; // atan2d not available until MATLAB 2012b

		// calculate altitude
		var alt = U * (1.0 - b2 / (a * V));

		var LLA = {lat: lat, lon: lon, alt: alt};

		return LLA;
	}
	
function ned_to_ecf(ned_value, orp_ecf){

	//var temp_ecf_value = create2DArray(3);

	var temp_ecf_value = transposeMatrix(ecf_ned_rot_mat(orp_ecf));

	var ecf_value = [];
	ecf_value = multiplyMatrix21(temp_ecf_value, ned_value);
	return ecf_value;
}




function transposeMatrix(orig_ecf_value) {

	var transposed_ecf_value = create2DArray(orig_ecf_value.length);
	// transpose
	if (orig_ecf_value.length > 0) {
		for (var i=0; i < orig_ecf_value.length; i++) { //assuming square matrix
			for (var j=0; j < orig_ecf_value.length; j++) {
				transposed_ecf_value[i][j] = orig_ecf_value[j][i];
			}
		}
	}
	return transposed_ecf_value;
}



function ecf_ned_rot_mat(orp_ecf) {

	var rad2deg = 180/Math.PI;
	var deg2rad = Math.PI/180;

	var orp_lla = ECEFtoLLA(orp_ecf[0], orp_ecf[1], orp_ecf[2]);

	var rot_mat1 = [[ Math.cos((-90-orp_lla.lat)*deg2rad), 0.0, -Math.sin((-90-orp_lla.lat)*deg2rad) ],
			[ 0.0, 1.0, 0.0 ],
			[ Math.sin((-90-orp_lla.lat)*deg2rad), 0.0, Math.cos((-90-orp_lla.lat)*deg2rad) ]];

	var rot_mat2 = [[ Math.cos((orp_lla.lon)*deg2rad), Math.sin((orp_lla.lon)*deg2rad), 0.0 ],
			[ -Math.sin((orp_lla.lon)*deg2rad), Math.cos((orp_lla.lon)*deg2rad), 0.0 ],
			[ 0.0, 0.0, 1.0 ]];

	var rot_mat = multiplyMatrix22(rot_mat1, rot_mat2);

	return rot_mat;
}

function multiplyMatrix21(first, second) {
	//2 dimensional matrix times a 1 dimensional matrix
	var m, n, p, q, c, d, k;
	var sum = 0;

	m = first.length; //num rows for first matrix;
	n = first[0].length; //num cols for first matrix;
	p = second.length; //num rows for second matrix;
	q = second.length; //num cols for second matrix;
	
	var multiplied = create2DArray(first.length);


	for ( c = 0 ; c < m ; c++ )
	{
		for ( d = 0 ; d < q ; d++ )
		{
			for ( k = 0 ; k < p ; k++ )
			{
				sum = sum + first[c][k]*second[k];
			}
			multiplied[c] = sum;
			sum = 0;
		}
	}
	return multiplied;
}

function multiplyMatrix22(first, second){
	//2 dimensional matrix times a 2 dimensional matrix
	var m, n, p, q, c, d, k;
	var sum = 0;

	m = first.length; //num rows for first matrix;
	n = first[0].length; //num cols for first matrix;
	p = second.length; //num rows for second matrix;
	q = second[0].length; //num cols for second matrix;

	var multiplied = create2DArray(first.length);

	for ( c = 0 ; c < m ; c++ )
	{
		for ( d = 0 ; d < q ; d++ )
		{
			for ( k = 0 ; k < p ; k++ )
			{
				sum = sum + first[c][k]*second[k][d];
			}
			multiplied[c][d] = sum;
			sum = 0;
		}
	}
	return multiplied;
}

function create2DArray(rows){
	var arr = [];
	for (var i=0;i<rows;i++){
		arr[i] = [];
	}
	return arr;
}
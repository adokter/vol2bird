/** Functions for reading NEXRAD WSR88-D data
 * @file librsl.c
 * @author Adriaan Dokter
 * @date 2017-05-05
 */

#ifdef RSL

#include "rsl.h"
#include <libgen.h>
#include "polarvolume.h"
#include "polarscan.h"
#include "constants.h"
#include "libvol2bird.h"
#include "rave_debug.h"

// non-public function prototypes (local to this file/translation unit)

PolarVolume_t* PolarVolume_vol2bird_RSL2Rave(Radar* radar, float rangeMax);

PolarVolume_t* PolarVolume_RSL2Rave(Radar* radar, float rangeMax);

PolarScan_t* PolarScan_RSL2Rave(Radar *radar, int iScan, float rangeMax);

PolarScanParam_t* PolarScanParam_RSL2Rave(Radar *radar, float elev, int RSL_INDEX,float rangeMax, double *scale);
    
int rslCopy2Rave(Sweep *rslSweep,PolarScanParam_t* scanparam);


// non-public function declarations (local to this file/translation unit)

// copies a RSL sweep to a Rave scan
int rslCopy2Rave(Sweep *rslSweep,PolarScanParam_t* scanparam){
    float value;
    double setvalue;
    float rscale;
    int rayindex=0;
    Ray *rslRay;
    long nrays,nbins;
    
    rslRay = RSL_get_first_ray_of_sweep(rslSweep);
    
    if (rslRay == NULL) return 0;
    
    nbins = PolarScanParam_getNbins(scanparam);
    nrays = PolarScanParam_getNrays(scanparam);

    if (nbins == 0 || nrays == 0) return 0;

    for(int iRay=0; iRay<rslSweep->h.nrays; iRay++){
        // determine at which ray index we are in the rave scanparam
        // adding half a ray bin width, to get into the middle of the ray bin
        rayindex=ROUND(nrays*(rslRay->h.azimuth+180.0/nrays)/360.0);
        // get the range gate size for this ray
        rscale = rslRay->h.gate_size;
        // only values between 0 and 360 degrees permitted
        if (rayindex >= nrays) rayindex-=nrays;
        // loop over range bins
        int iBinStart = ROUND((rslRay->h.range_bin1 + 0.5*rscale)/rscale);
        for(int iBin=iBinStart; iBin<nbins; iBin++){
            value=RSL_get_value_from_ray(rslRay, iBin*rscale/1000);
            if (value == BADVAL || value == RFVAL){
                // BADVAL is used in RSL library to encode for undetects, but also for nodata
                // In most cases we are dealing with undetects, so encode as such
                // RFVAL is a range folded value, see RSL documentation.
                setvalue=PolarScanParam_getUndetect(scanparam);
            }
            else{
                setvalue=(value-PolarScanParam_getOffset(scanparam))/PolarScanParam_getGain(scanparam);
            }
            PolarScanParam_setValue(scanparam, iBin, rayindex, setvalue);
        }
        rslRay=RSL_get_next_cwise_ray(rslSweep, rslRay);
    }
    
    return 1;
}


PolarScanParam_t* PolarScanParam_RSL2Rave(Radar *radar, float elev, int RSL_INDEX,float rangeMax, double *scale){
    Volume *rslVolume;
    Sweep *rslSweep;
    Ray *rslRay;
    PolarScanParam_t* param = NULL;
    
    float offset,gain;
    int nbins,nrays,rscale;
    
    char* DBZH="DBZH";
    char* VRADH="VRADH";
    char* RHOHV="RHOHV";
    char* WRADH="WRADH";
    char* TH="TH";
    char* ZDR="ZDR";
    char* PHIDP="PHIDP";
    char* KDP="KDP";
    char* VRAD2="VRAD2";
    char* VRAD3="VRAD3";
    char *name;

    if(radar == NULL) {
        fprintf(stderr, "Warning: RSL radar object is empty...\n");
        return param;
    }
    
    rslVolume = radar->v[RSL_INDEX];

    if(rslVolume == NULL) {
        fprintf(stderr, "Warning: RSL volume %i not found by PolarScanParam_RSL2Rave...\n",RSL_INDEX);
        return param;
    }

    rslSweep = RSL_get_sweep(rslVolume, elev);

    if(rslSweep == NULL) {
        fprintf(stderr, "Warning: RSL sweep of volume %i not found by PolarScanParam_RSL2Rave...\n",RSL_INDEX);
        return param;
    }

    rslRay = RSL_get_first_ray_of_sweep(rslSweep);
    

    if(rslRay == NULL) {
        fprintf(stderr, "Warning: RSL first ray of volume %i not found by PolarScanParam_RSL2Rave...\n",RSL_INDEX);
        return param;
    }
    
    if(ABS(rslSweep->h.elev-elev)>ELEVTOL){
        fprintf(stderr,"Warning: elevation angle mistmatch in PolarScanParam_RSL2Rave (requested %f, found %f)...\n",elev,rslSweep->h.elev);
        return param;
    }
    
    switch (RSL_INDEX) {
        case DZ_INDEX : 
            name = DBZH;
            offset = 0;
            gain = 1;
            break;
        case VR_INDEX :
            name = VRADH;
            offset = 0;
            gain = 1;
            break;        
        case RH_INDEX :
            name = RHOHV;
            offset = 0;
            gain = 1;
            break;
        case SW_INDEX :
            name = WRADH;
            offset = 0;
            gain = 1;
            break;
        case ZT_INDEX :
            name = TH;
            offset = 0;
            gain = 1;
            break;
        case DR_INDEX :
            name = ZDR;
            offset = 0;
            gain = 1;
            break;
        case PH_INDEX :
            name = PHIDP;
            offset = 0;
            gain = 1;
            break;
        case KD_INDEX :
            name = KDP;
            offset = 0;
            gain = 1;
            break;
        case V2_INDEX :
            name = VRAD2;
            offset = 0;
            gain = 1;
            break;
        case V3_INDEX :
            name = VRAD3;
            offset = 0;
            gain = 1;
            break;
        default :
            fprintf(stderr, "Something went wrong; RSL scan parameter not implemented in PolarScanParam_RSL2Rave\n");
            return param;

    }
 
    // get the number of range bins
    rscale = rslRay->h.gate_size;
    nbins = 1+rslRay->h.nbins+rslRay->h.range_bin1/rscale;
    nbins = MIN(nbins, ROUND(rangeMax/rscale));
    // return rscale through the scale parameter
    *scale = (double) rscale;
        
    // estimate the number of azimuth bins
    // early WSR88D scans have somewhat variable azimuth spacing
    // therefore we resample the data onto a regular azimuth grid
    // azimuth bin size is round off to 1/n with n positive integer
    // i.e. either 1, 0.5, 0.25 degrees etc.
    nrays = 360*ROUND((float) rslSweep->h.nrays/360.0);
    if (nrays != rslSweep->h.nrays){
        fprintf(stderr, "Warning: resampling %s sweep at elevation %f (%i rays into %i azimuth-bins) ...\n",name,elev,rslSweep->h.nrays,nrays);
    }

    param = RAVE_OBJECT_NEW(&PolarScanParam_TYPE);
    PolarScanParam_setQuantity(param, name);
    PolarScanParam_createData(param,nbins,nrays,RaveDataType_DOUBLE);
    PolarScanParam_setOffset(param,offset);
    PolarScanParam_setGain(param,gain);
    PolarScanParam_setNodata(param,RSL_NODATA);
    PolarScanParam_setUndetect(param,RSL_UNDETECT);

    // initialize the data field
    for(int iRay=0; iRay<nrays; iRay++){
        for(int iBin=0; iBin<nbins; iBin++){
            PolarScanParam_setValue(param, iBin, iRay, PolarScanParam_getNodata(param));
        }
    }

    // Fill the PolarScanParam_t objects with corresponding RSL data
    //XXXX Check that the encoding works correctly
    rslCopy2Rave(rslSweep,param);
    
    return param;
}


PolarScan_t* PolarScan_RSL2Rave(Radar *radar, int iScan, float rangeMax){
    
    PolarScanParam_t* param;
    PolarScan_t* scan = NULL;
    Volume *rslVol = NULL;
    Ray *rslRay = NULL;
    double elev = 0;
    double nyq_vel = 0;
    double rscale;
    
    if(radar == NULL) {
        fprintf(stderr, "Error: RSL radar object is empty...\n");
        return scan;
    }
    
    // get first RSL volume
    for (int iParam = 0; iParam < radar->h.nvolumes; iParam++){
        if(radar->v[iParam] == NULL) continue;
        rslVol = radar->v[iParam];    
        elev = rslVol->sweep[iScan]->h.elev;    
        break;
    }
    
    if(rslVol == NULL) {
        fprintf(stderr, "Error: RSL radar object is empty...\n");
        return scan;
    }

    if(iScan>rslVol->h.nsweeps-1) {
        fprintf(stderr, "Error: iScan larger than # sweeps...\n");
        return scan;
    }

    // all checks on sweep passed
    // start a new scan object
    scan = RAVE_OBJECT_NEW(&PolarScan_TYPE);
                
    // add attribute Elevation to scan
    PolarScan_setElangle(scan, (double) rslVol->sweep[iScan]->h.elev*PI/180);

    // add attribute Beamwidth to scan
    PolarScan_setBeamwidth(scan, (double) rslVol->sweep[iScan]->h.beam_width);

    // add attribute Nyquist velocity to scan
    rslRay = RSL_get_first_ray_of_volume(rslVol);
    nyq_vel = rslRay->h.nyq_vel;
    
    // if no nyquist velocity found, try it with the native RSL function
    if(nyq_vel == 0){
        nyq_vel = RSL_get_nyquist_from_radar(radar);
    }
    
    RaveAttribute_t* attr_NI = RaveAttributeHelp_createDouble("how/NI", (double) nyq_vel);
    if (attr_NI == NULL || nyq_vel == 0){
        fprintf(stderr, "warning: no valid Nyquist velocity found in RSL polar volume\n");
    }
    else{    
        PolarScan_addAttribute(scan, attr_NI);
    }

    // add range scale Atribute to scan
    rscale = rslRay->h.gate_size;
    PolarScan_setRscale(scan, rscale);

    // loop through the volume pointers
    // iParam gives you the XX_INDEX flag, i.e. scan parameter type
    int result = 0;
    double scale = 0;
    for (int iParam = 0; iParam < radar->h.nvolumes; iParam++){
        
        if(radar->v[iParam] == NULL) continue;
        
        param = PolarScanParam_RSL2Rave(radar, elev, iParam, rangeMax, &scale);
        if(param == NULL){
            fprintf(stderr, "PolarScanParam_RSL2Rave returned empty object for parameter %i\n",iParam); 
            goto done;
        }
        
        // add the scan parameter to the scan, assuming their dimensions fit
        // for NEXRAD legacy data this will not work, as dimensions differ
        result = PolarScan_addParameter(scan, param);
        
        if(result == 0){
            fprintf(stderr, "Warning: dimensions of scan parameter %i do not match scan dimensions, resampling ...\n",iParam);
            
            PolarScanParam_t *param_proj;
            // project the scan parameter on the grid of the scan
            param_proj = PolarScanParam_project_on_scan(param, scan, scale);
            // try to add it again
            result = PolarScan_addParameter(scan, param_proj);
            
            if(result == 0){
                fprintf(stderr, "PolarScan_RSL2Rave failed to add parameter %i to RAVE polar scan\n",iParam); 
                RAVE_OBJECT_RELEASE(param_proj);
            }
 
            RAVE_OBJECT_RELEASE(param);
        }
        
    }

    done:
        return scan;
}

// maps a RSL polar volume to a RAVE polar volume NEW NEW NEW
PolarVolume_t* PolarVolume_RSL2Rave(Radar* radar, float rangeMax){
        
    // the RAVE polar volume to be returned by this function
    PolarVolume_t* volume = NULL;
    
    if(radar == NULL) {
        fprintf(stderr, "Error: RSL radar object is empty...\n");
        return volume;
    }

    // sort the scans (sweeps) and rays
    if(RSL_sort_radar(radar) == NULL) {
        fprintf(stderr, "Error: failed to sort RSL radar object...\n");
        goto done;
    }
    
    Volume *rslVol = NULL;
    Ray* rslRay = NULL;
    PolarScan_t* scan = NULL;

    // several checks that the volumes contain data -- should be outside this function!
    // should have a specific vol2bird version, and a rsl2odim version
    // this should be in vol2bird version only

    // find the first non-zero RSL volume
    for (int iParam = 0; iParam < radar->h.nvolumes; iParam++){
        if(radar->v[iParam] == NULL) continue;
        rslVol = radar->v[iParam];        
        break;
    }

    // find the largest shared maximum range
    float maxRange=FLT_MAX;
    float iRange;
    for (int iParam = 0; iParam < radar->h.nvolumes; iParam++){
        if(radar->v[iParam] == NULL) continue;
        rslRay = RSL_get_first_ray_of_volume(radar->v[iParam]);
        iRange=rslRay->h.range_bin1 + rslRay->h.nbins * rslRay->h.gate_size;
        if(iRange < maxRange) maxRange=iRange;
    }
    // if largest shared maximum range is larger than requested rangeMax, use the requested value
    if(rangeMax<maxRange) maxRange=rangeMax;

    // retrieve the first ray, for extracting some metadata
    rslRay = RSL_get_first_ray_of_volume(rslVol);
    if (rslRay == NULL){
        fprintf(stderr, "Error: RSL radar object contains no rays...\n");
        goto done;
    }
        
    // all checks on RSL object passed
    // make a new rave polar volume object
    volume = RAVE_OBJECT_NEW(&PolarVolume_TYPE);
    
    if (volume == NULL) {
        RAVE_CRITICAL0("Error: failed to create polarvolume instance");
        goto done;
    }

    // add attribute data to RAVE polar volume
    // first, copy metadata stored in radar header
    char pvtime[7];
    char pvdate[9];
    char *pvsource = malloc(strlen(radar->h.name)+strlen(radar->h.city)+strlen(radar->h.state)+strlen(radar->h.radar_name)+30);
    sprintf(pvtime, "%02i%02i%02i",radar->h.hour,radar->h.minute,ROUND(radar->h.sec));
    sprintf(pvdate, "%04i%02i%02i",radar->h.year,radar->h.month,radar->h.day);
    sprintf(pvsource, "RAD:%s,PLC:%s,state:%s,radar_name:%s",radar->h.name,radar->h.city,radar->h.state,radar->h.radar_name);
    fprintf(stderr,"Reading RSL polar volume with nominal time %s-%s, source: %s\n",pvdate,pvtime,pvsource);    
    PolarVolume_setTime(volume,pvtime);
    PolarVolume_setDate(volume,pvdate);
    PolarVolume_setSource(volume,pvsource);
    PolarVolume_setLongitude(volume,(double) (radar->h.lond + radar->h.lonm/60.0 + radar->h.lons/3600.0)*PI/180);
    PolarVolume_setLatitude(volume,(double) (radar->h.latd + radar->h.latm/60.0 + radar->h.lats/3600.0)*PI/180);
    PolarVolume_setHeight(volume, (double) radar->h.height);

    // second, copy volume coverage pattern (VCP) information to how/vcp attribute (NEXRAD specific)
    int vcp = radar->h.vcp;
    RaveAttribute_t* attr_vcp = RaveAttributeHelp_createLong("how/vcp", (long) vcp);
    if (attr_vcp == NULL){
        fprintf(stderr, "warning: no valid VCP value found in RSL polar volume\n");
    }
    else{    
        PolarVolume_addAttribute(volume, attr_vcp);
    }

    // third, copy metadata stored in ray header; assume attributes of first ray applies to entire volume
    float wavelength = rslRay->h.wavelength*100;
    RaveAttribute_t* attr_wavelength = RaveAttributeHelp_createDouble("how/wavelength", (double) wavelength);
    if (attr_wavelength == NULL && wavelength > 0){
        fprintf(stderr, "warning: no valid wavelength found in RSL polar volume\n");
    }
    else{    
        PolarVolume_addAttribute(volume, attr_wavelength);
    }
    
    // read the RSL scans (sweeps) and add them to RAVE polar volume
    int result;
    for (int iScan = 0; iScan < rslVol->h.nsweeps; iScan++){
        scan = PolarScan_RSL2Rave(radar, iScan, maxRange);
        // Add the scan to the volume
        result = PolarVolume_addScan(volume,scan);
        if(result == 0){
           fprintf(stderr, "PolarVolume_RSL2Rave failed to add RSL scan %i to RAVE polar volume\n",iScan); 
        }
    }
   
    free(pvsource);
    
    done:
        return volume;
}


// has extra checks in place ... XXX
PolarVolume_t* PolarVolume_vol2bird_RSL2Rave(Radar* radar, float rangeMax){
    // pointers to RSL polar volumes of reflectivity, velocity, respectively.
    Volume *rslVolZ,*rslVolV;
    // pointer to a RSL ray

    // retrieve the polar volumes from the Radar object
    rslVolZ = radar->v[DZ_INDEX];
    rslVolV = radar->v[VR_INDEX];
    
    PolarVolume_t* volume = NULL;
    
    // several checks that the volumes contain data
    if (rslVolZ == NULL){
        fprintf(stderr, "Error: RSL radar object contains no reflectivity volume...\n");
        goto done;
    }
    else if (rslVolV == NULL){
        fprintf(stderr, "Error: RSL radar object contains no radial velocity volume...\n");
        goto done;
    }
    
    volume = PolarVolume_RSL2Rave(radar, rangeMax);
    
    done:
        return volume;
}


PolarVolume_t* vol2birdGetRSLVolume(char* filename, float rangeMax, int small) {
    Radar *radar;
    PolarVolume_t* volume = NULL;

    // if small, only read reflectivity, velocity, Rho_HV        
    // else select all scans
    if(small) RSL_select_fields("dz","vr","sw","rh", NULL);
    else RSL_select_fields("dz","vr","sw","zt","dr","rh","ph","kd", NULL);
    
    RSL_read_these_sweeps("all",NULL);
    
    // read the file to a RSL radar object
    
    // according to documentation of RSL it is not required to parse a callid
    // but in practice it is for WSR88D.
    char* base = basename(filename);
    char callid[5];
    strncpy(callid, base,4);
    callid[4] = 0; //null terminate destination
    radar = RSL_anyformat_to_radar(filename,callid);

    if (radar == NULL) {
        fprintf(stderr, "critical error, cannot open file %s\n", filename);
        return NULL;
    }
    
    // convert RSL object to RAVE polar volume
    
    volume = PolarVolume_vol2bird_RSL2Rave(radar, rangeMax);
    
    RSL_free_radar(radar);
    
    return(volume);
    
}

#endif

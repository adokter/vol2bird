import unittest
import math
import numpy as np
import os
import _raveio
import _verticalprofile
import rave_pgf_vol2bird_plugin

class PyVol2BirdTest(unittest.TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def _dim(self, a):
        if isinstance(a, list):
            return [len(a)] + self._dim(a[0])
        return []


    def _assert_attributes_covered_by(self, refdsc, ref, testdsc, test):
        '''Call with either VerticalProfileCore objects or RaveFieldCore objects'''
        for attrname in ref.getAttributeNames():
            refval = ref.getAttribute(attrname)
            testval = test.getAttribute(attrname)
            self.assertIsNotNone(testval, "%s is missing attribute %s" % (testdsc, attrname))
            if attrname not in ["what/nodata", "what/undetect"]:
                self.assertEqual(refval, testval, "Attribute %s values mismatch between %s and %s"
                    % (attrname, refdsc, testdsc))


    def _assert_is_covered_by(self, refdsc, refvpr, testdsc, testvpr):
        self._assert_attributes_covered_by(refdsc, refvpr, testdsc, testvpr)
        for reffield in refvpr.getFields():
            reffieldname = reffield.getAttribute("what/quantity")
            testfield = testvpr.getField(reffieldname)
            self.assertIsNotNone(testfield, "%s does not have field %s" % (testdsc, reffieldname))
            self._assert_attributes_covered_by(refdsc, reffield, testdsc, testfield)


    def _data_matches(self, field, attrname):
        '''Return a bool array with True wherever the field data matches the attribute.'''
        retval = None
        data = field.getData()
        if (attrname in field.getAttributeNames()):
            value = field.getAttribute(attrname)

            if math.isnan(value):
                # compare NaN as if it were a normal value
                retval = np.isnan(data)
            else:
                retval = (data == value)

        else:
            # attribute does not exist, so no values match it
            # create an array of Falses of same shape as field
            retval = np.zeros(data.shape, dtype=np.bool_)

        return retval


    def _implies(self, test, ref):
        '''Test says it's true, ref says it's false, conflict!

           We use the "implies" logical operator to compare the test and reference, as per
           the following table.

           Test -> Reference

           False -> False returns True
           False -> True  returns True
           True  -> False returns False
           True  -> True  returns True

           As a result, a mismatch is detected if the test field contains e.g. a nodata
           value, while the reference field contains a value that can not be interpreted
           as nodata.

           If the reference file has what/nodata set to NaN, and also what/undetect set
           to NaN, and a NaN in position [0, 0], then a test file with a nodata value at
           [0, 0] will be accepted, as will a test file with an undetect value at [0, 0].

           Conversely, if the test file has what/nodata set to NaN, and also what/undetect
           set to NaN, and a NaN in position [0, 0], then it will only be accepted if the
           reference file also has what/nodata and what/undetect set to the same value,
           and has that value at position [0, 0].
        '''
        return np.logical_and(test, np.logical_not(ref))


    def _any_of(self, nparray1, nparray2, nparray3):
        return np.logical_or(
            np.logical_or(nparray1, nparray2),
            nparray3)


    def _assert_contents_equivalent(self, refdsc, refvpr, testdsc, testvpr):
        for reffield in refvpr.getFields():
            # look up matching field in vertical profile under test
            reffieldname = reffield.getAttribute("what/quantity")
            testfield = testvpr.getField(reffieldname)

            # get field data
            refdata = reffield.getData();
            testdata = testfield.getData();

            # check that fields have the same shape
            self.assertEqual(refdata.shape, testdata.shape,
                "Shapes of field %s mismatch between %s and %s" %
                    (reffieldname, testdsc, refdsc))

            # categorize values in reference field
            ref_isnodata = self._data_matches(reffield, "what/nodata")
            ref_isundetect = self._data_matches(reffield, "what/undetect")
            ref_isundefined = self._data_matches(reffield, "what/undefined")
            ref_isdata = np.logical_not(self._any_of(ref_isnodata, ref_isundetect, ref_isundefined))

            # categorize values in field under test
            test_isnodata = self._data_matches(testfield, "what/nodata")
            test_isundetect = self._data_matches(testfield, "what/undetect")
            test_isundefined = self._data_matches(testfield, "what/undefined")
            test_isdata = np.logical_not(self._any_of(test_isnodata, test_isundetect, test_isundefined))

            # check categories against each other
            isnodata_mismatch = self._implies(test_isnodata, ref_isnodata)
            isundetect_mismatch = self._implies(test_isundetect, ref_isundetect)
            isundefined_mismatch = self._implies(test_isundefined, ref_isundefined)
            isdata_mismatch = self._implies(test_isdata, ref_isdata)

            # check values wherever we have data
            value_mismatch = (np.logical_and(ref_isdata, (refdata != testdata))).any()

            # ensure everything is okay
            msg = ("%s has a %%(attr)s in field %s where %s says it cannot be %%(attr)s" %
                (testdsc, reffieldname, refdsc))
            self.assertFalse(isnodata_mismatch.any(), msg % {"attr": "nodata"})
            self.assertFalse(isundetect_mismatch.any(), msg % {"attr": "undetect"})
            self.assertFalse(isundefined_mismatch.any(), msg % {"attr": "undefined"})
            self.assertFalse(isdata_mismatch.any(), msg % {"attr": "data"})
            self.assertFalse(value_mismatch.any(), "Values of %s and %s differ in field %s" %
                (testdsc, refdsc, reffieldname))


    def _assert_equal(self, refdsc, refvpr, testdsc, testvpr):
        self._assert_is_covered_by(refdsc, refvpr, testdsc, testvpr)
        self._assert_is_covered_by(testdsc, testvpr, refdsc, refvpr)
        self._assert_contents_equivalent(refdsc, refvpr, testdsc, testvpr)


    def _loadProfile(self, filename):
        file = _raveio.open(filename)
        object = file.object
        self.assertTrue(_verticalprofile.isVerticalProfile(object))
        return object


    def test_Vol2Bird(self):
        files = os.listdir(os.path.join("fixtures", "vol2bird", "pvol"))
        for testname in files:
            testfilename = os.path.join("fixtures", "vol2bird", "pvol", testname)
            v2boutname = rave_pgf_vol2bird_plugin.generate([testfilename], [])

            testvp = self._loadProfile(v2boutname)
            testdsc = "vol2bird output of %s" % testname

            refname = testname.replace("_pvol_", "_vp_")
            reffilename = os.path.join("fixtures", "vol2bird", "vp", refname)
            refvp = self._loadProfile(reffilename)

            self._assert_equal(refname, refvp, testdsc, testvp)
            os.remove(v2boutname)



** Authors: David Aikman and Simone Maso 
** Date: May 2025

graph close _all
clear
cls

* set the directory
global	base_dir "C:\Users\k2370179\Dropbox\BoE-KCL Macro Forecasting\Data"

*=========================*
** from daily to monthly **
*=========================*

* import the series
import excel using "${base_dir}/others/countercyclical-capital-buffer.xlsx", sheet("7. Global corp. bond spreads") cellrange(A4:E7372) firstrow clear

* clean 
rename *, lower 
keep date gbpinvestmentgrade

* generate monthly date
gen mdate = mofd(date)
format mdate %tm

* compute the average 
bys mdate: egen avg_corp_spread = mean(gbpinvestmentgrade)

* final clean 
keep mdate avg_corp_spread
bys *: keep if _n == 1

*=============================*
** from monthly to quarterly **
*=============================*

*************************************
* BOND SPREAD AMBROGIO-CESA-BIANCHI *
*************************************

* import the series
import excel using "${base_dir}/others/ambrogio_cesa_bianchi_data_series.xlsx", sheet("UK") cellrange(A1:B487) firstrow clear

* clean the data variable 
gen date = qofd(dofm(monthly(month, "YM")))
format date %tq

* gen quarterly average 
bys date: egen bond_spread_q = mean(CS_i_10YR)

* clean
keep date bond_spread_q
bys *: keep if _n == 1

* copy paste in GaRDataRaw.xlsx

*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

*******************************************
* BOND SPREAD "A MILLENIUM OF MACRO DATA" *
*******************************************

* import the series
import excel using "${base_dir}/IaRDataRaw_monthly.xlsx", sheet("bond_spread") cellrange(A1:B666) firstrow clear

* clean the data variable 
gen date = qofd(dofm(monthly(A, "YM")))
format date %tq

* gen quarterly average 
bys date: egen bond_spread_q = mean(UK)

* clean
keep date bond_spread_q
bys *: keep if _n == 1

* copy paste in GaRDataRaw_quarterly.xlsx

*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

********
* PMIs *
********

* conver PMI 

* import the series
import excel using "${base_dir}/Alternative_variables2.xlsx", sheet("PMI") cellrange(A1:C427) firstrow clear

* clean the data variable 
gen date = qofd(Date)
format date %tq

foreach var in cipsto cipste {
	* gen quarterly average 
	bys date: egen `var'_q = mean(`var')
}

* clean
keep date *_q
bys *: keep if _n == 1

* copy paste in GaRDataRaw_quarterly.xlsx

*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

**********************
* CONSENSUS FORECAST *
**********************

* import the series
import excel using "${base_dir}/Alternative_variables2.xlsx", sheet("CONSENSUS") cellrange(A1:C400) firstrow clear

* clean the data variable 
gen date = qofd(Date)
format date %tq

foreach var in rgdp0_m	rgdp1_m	 {
	* gen quarterly average 
	bys date: egen `var'_q = mean(`var')
}

* clean
keep date *_q
bys *: keep if _n == 1

* copy paste in GaRDataRaw_quarterly.xlsx

*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

*********************
* YIELD CURVE SLOPE *
*********************

* Create a temp file to store the master dataset
tempfile master

* Track first iteration
local first = 1

* Import the data
foreach yy in 1979 1985 1990 1995 2000 2005 2016 2025 {

    * Define the corresponding end year
    local enddate = cond(`yy' == 1979, "1984", ///
                    cond(`yy' == 1985, "1989", ///
                    cond(`yy' == 1990, "1994", ///
                    cond(`yy' == 1995, "1999", ///
                    cond(`yy' == 2000, "2004", ///
                    cond(`yy' == 2005, "2015", ///
                    cond(`yy' == 2016, "2024", ///
                    cond(`yy' == 2025, "present", ""))))))))

    * Set the correct cell range based on the year
    local cellrange = cond(`yy' == 1979, "A5:AY1571", ///
                    cond(`yy' == 1985, "A5:AY1309", ///
                    cond(`yy' == 1990, "A5:AY1310", ///
                    cond(`yy' == 1995, "A5:AY1310", ///
                    cond(`yy' == 2000, "A5:AY1310", ///
                    cond(`yy' == 2005, "A5:AY2874", ///
                    cond(`yy' == 2016, "A5:CC2353", ///
                    cond(`yy' == 2025, "A5:CC157", ""))))))))

    * Set the correct sheet name
    local sheetname = cond(inlist(`yy', 2005, 2016, 2025), "4. spot curve", "4. nominal spot curve")

    * Import the Excel file with the correct range and sheet
    import excel using "${base_dir}/others/yield_curve_data/GLC Nominal daily data_`yy' to `enddate'.xlsx", ///
        sheet("`sheetname'") cellrange(`"`cellrange'"') firstrow clear

		
   * Save or append
    if `first' {
        save `master'
        local first = 0
    }
    else {
        append using `master'
        save `master', replace
    }
}

* sort the data 
sort date

* keep 1 year and 10 year 
keep date mat1 mat10
gen slope = mat10 - mat1

* generate average in the quarter 
gen qdate = qofd(date)
format qdate %tq
bys qdate: egen avg_slope = mean(slope)

* remove dups 
keep qdate avg_slope
bys *: keep if _n == 1

* copy paste in GaRDataRaw_quarterly.xlsx


*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

**************************
* INFLATION EXPECTATIONS *
**************************

* from monthly to quarterly 
import excel using "${base_dir}/IaRDataRaw_monthly.xlsx", sheet("infl1_m") cellrange(A1:B666) firstrow clear

* clean the data variable 
gen date = qofd(dofm(monthly(A, "YM")))
format date %tq

* gen quarterly average 
bys date: egen inf_exp_q = mean(UK)

* clean
keep date inf_exp_q
bys *: keep if _n == 1

* copy paste in GaRDataRaw_quarterly.xlsx

*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

***************
* MONTHLY GDP *
***************

* import the data Gross Value Added - Monthly (Index 1dp) :CVM SA ECY2
import excel using "${base_dir}/others/mgdp.xlsx", sheet("data") cellrange(A7:B348) firstrow clear

* compute the year-on-year growth between the the average of the first two months of the current quarter and the average of the 3 months in the same quarter the previous year 
gen mdate = monthly(ImportantNotes, "YM")
format mdate %tm
gen month_num = month(dofm(monthly(ImportantNotes, "YM")))

* rename 
drop ImportantNotes
rename B mgdp

* generate mdate with only the first two monhts in the quarter
gen mgdp2m = mgdp
replace mgdp2m = . if inlist(month, 3, 6, 9, 12)

* get the average in the quarter for mgdp2m and mgdp
gen qdate = qofd(dofm(mdate))
format qdate %tq
collapse (mean) mgdp mgdp2m, by(qdate)

* get the y-o-y growth 
sort qdate 
gen growth = mgdp2m/mgdp[_n-4] - 1

* copy paste in GaRDataRaw_quarterly.xlsx


*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

**************
* EPU INDEX  *
**************

* import the data
import excel using "${base_dir}/others/EPU_index.xlsx", sheet("data") cellrange(A1:F1508) firstrow clear

* clean the data variable 
gen qdate = qofd(date)
format qdate %tq

* keep date and splicedEPU 
keep qdate splicedEPU

* compute quarterly averages
bys qdate: egen avg_EPU = mean(splicedEPU)
keep qdate avg_EPU
bys *: keep if _n == 1

* copy paste in GaRDataRaw_quarterly.xlsx

*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

*******************************
* LPMVTVX - housing approvals *
*******************************

* import the data
import excel using "${base_dir}/others/housing_approval.xlsx", sheet("Sheet1") cellrange(A2:C383) firstrow clear

* clean the data variable 
gen qdate = qofd(A)
format qdate %tq

* keep date and housign appr  
keep qdate Monthlynumberoftotalsterling

* compute quarterly averages
bys qdate: egen housing_app = sum(Monthlynumberoftotalsterling)
keep qdate housing_app
bys *: keep if _n == 1

* copy paste in GaRDataRaw_quarterly.xlsx

*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

********************
* IMPORT DEFLATOR  *
********************

* linearly interpolate the month within a quarter

* import the data
import excel using "${base_dir}/others/import_deflator.xlsx", sheet("import_deflator") cellrange(A86:B366) clear

* rename the col  
rename (A B) (date qseries)

* extend the series
gen qdate = quarterly(date, "YQ")
format qdate %tq
gen mdate = mofd(dofq(qdate)) + 2   // assign the value of the quarter to the last month in the quarter
format mdate %tm

* Declare time series structure
tsset mdate, monthly

* Expand: fill in all monthly slots (adds 2 missing months per quarter)
tsfill

* clean the environment 
drop qdate date 

* 4. Interpolate linearly for missing months
ipolate qseries mdate, gen(mseries)

* copy paste in GaRDataRaw_quarterly.xlsx

*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

% GETDBDATA: An automatically-generated MATLAB script for reading
% historical credit rating data from a database.

% Set preferences with setdbprefs.
s.DataReturnFormat = 'structure';
s.ErrorHandling = 'store';
s.NullNumberRead = 'NaN';
s.NullNumberWrite = 'NaN';
s.NullStringRead = 'null';
s.NullStringWrite = 'null';
s.JDBCDataSourceFile = '';
s.UseRegistryForSources = 'yes';
s.TempDirForRegistryOutput = 'C:\Temp';
s.DefaultRowPreFetch = '10000';
setdbprefs(s)

% Make connection to database.  Note that the password has been omitted.
% Using ODBC driver.
conn = database('Historical Credit Ratings','','password');

% Read data from database.
e = exec(conn,'SELECT ALL ID,WC_TA,RE_TA,EBIT_TA,MVE_BVTD,S_TA,Industry,Rating FROM Ratings');
e = fetch(e);
close(e)

% Close database connection.
close(conn)

% Create a MATLAB structure from the returned data.
historicalData = e.data;

clear conn e s
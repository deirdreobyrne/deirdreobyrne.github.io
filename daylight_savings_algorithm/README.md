# Calculating Daylight Savings Time changeovers

## DST rules

There are [a limited number of rules](https://en.wikipedia.org/wiki/Daylight_saving_time_by_country) (16 at the time of writing) for when daylight savings
time comes into and out of effect. With the exception of the rules for Iran, the rules can be summarised as `[on|Fri before] [1st|2nd|last] [Thu|Fri|Sat|Sun]
in [Mar|Apr|Oct|Nov] at [time]`. Note, however, that DST rules can change, so it is prudent for a programmer to extend the possibilities, within reason. It
is also prudent to account for the fact that not all DST changes are plus one hour - at present, Lord Howe Island has a 30 minute DST change. DST has also
historically been 2 hours, especially during wartime.

## Step 1 - an algorithm to convert a date to a time interval since 1970

This algorithm is based on chapter 7 of Jean Meeus’ excellent 1991 book “Astronomical Algorithms” published by Willmann-Bell of Richmond, Virginia, USA.

```
// y - the year (>1970)
// m - the month (0 .. 11)
// d - the day of the month (1 .. 31)
// returns the number of days elapsed since 1970 Jan 1
int dayNumber(int y, int m, int d) {
  int days;
  if (m < 2) {
    y--;
    m+=12;
  }
  days = (int)(y/100);
  days = 365*y + (y>>2) - days + (days>>2) + 30*m + (int)((3*m+6)/5) + d - 719531;
  return days;
}
```
The algorithm works by first re-framing the start of the year as 1st March. This makes it easy to deal with leap years, as a leap year just "starts" one day
later than normal. Then we calculate the number of the century, and store the result temporarily in the `days` variable. Finally we calculate the actual number of
days elapsed.

The elements of that final formula are -
- `365*y + (y>>2)` - first step in calculating the number of days - assuming every 4th year is a leap year
- `- days + (days>>2)` - remembering that `days` here actually refers to the century number, this calculates the Gregorian correction to the calendar
- `30*m + (int)((3*m+6)/5)` - calculate the number of days elapsed to the start of the month. This formula neatly accounts for the varying lengths of
the months when `m` is 2 for March and 13 for February. Note that the `((3*m+6)/5)` needs to be calculated using integer arithmetic. Meeus actually gives
this part of the formula as `floor(30.6001*(m+2))`, which is safer if your computer insists on floating-point.
- `d` - add the day of the month
- ` - 719531` - a correction so that 1970 Jan 1 is zero

## Step 2 - Algorithm to implement the rules determining when the change occurs

Most computer languages have an implementation of something like `Date()`, from which the current year can be determined (e.g. by calling `Date().getYear()` or
equivalent), and the current number of milliseconds since 1970 Jan 1 (e.g. by calling `Date().getMillis()` or equivalent).

```
// year - the year for which we are determining the start/end of DST
// dstDaysPrior - if the DST rule is "Fri before last Sun", then this is 2, otherwise it is zero
// dstDayNumber - 0 for 1st, 1 for 2nd, and 4 for last
// dstDayOfWeek - 0 for Sunday
// dstMonth - 0 for January
// dstTime - the time of day (in hours) that the changeover happens
// dstTimezone - the "normal" timezone of the user in hours (positive east)
// dstOffset - the number of hours the DST changeover involves
// isStart - a boolean indicating whether we are looking at the start (as opposed to end) of DST
// Returns the number of milliseconds since 1970 that the change happens
int getDstChangeTime(year, dstDaysPrior, dstDayNumber, dstDayOfWeek, dstMonth, dstTime, dstTimezone, dstOffset, isStart) {
  int result;
  if (dstDayNumber == 4) { // last X of this month? Work backwards from 1st of next month.
    if (++dstMonth > 11) {
      year++;
      dstMonth-=12;
    }
  }
  result = dayNumber(year, dstMonth, 1); // 1970 Jan 1 was Thursday, so (result % 7) is 0 for Thursday, hence ((result + 4) % 7) is 0 for Sunday
  if (dstDayNumber == 4) {
    result -= 7 - (7 - ((result + 4) % 7) + dstDayOfWeek) % 7;
  } else {
    result += 7 * dstDayNumber + (14 + dstDayOfWeek - ((result + 4) % 7)) % 7;
  }
  result -= dstDaysPrior;
  result = (result * 86400) + (dstTime - dstTimezone - (isStart ? 0 : dstOffset)) * 3600;
  return result * 1000;
}
```






```
int getNextDSTChange(date, settings) {
  if (settings.has_dst) {
    var start = dstChangeTime(now.getFullYear(), settings.tz, settings.dst_start);
    var end = dstChangeTime(now.getFullYear(), settings.tz + settings.dst_size, settings.dst_end);
    if (start <= now.getTime()) {
      if (end <= now.getTime()) {
        // Both changes have happened for this year
        if (start < end) {
          // The start of DST is earlier than the end, so next change is a start of DST
          next_dst_change = { millis: dstChangeTime(now.getFullYear()+1, settings.tz, settings.dst_start), offset: settings.tz + settings.dst_size, is_start: true };
          setEffectiveTimezone(settings.tz);
        } else {
          // The end of DST is earlier than the start, so the next change is an end of DST
          next_dst_change = { millis: dstChangeTime(now.getFullYear()+1, settings.tz + settings.dst_size, settings.dst_end), offset: settings.tz, is_start: false };
          setEffectiveTimezone(settings.tz + settings.dst_size);
        }
      } else {
        next_dst_change = { millis: end, offset: settings.tz, is_start: false };
        setEffectiveTimezone(settings.tz + settings.dst_size);
      }
    } else {
      next_dst_change = { millis: start, offset: settings.tz + settings.dst_size, is_start: true };
      setEffectiveTimezone(settings.tz);
    }
    next_dst_change.show_icon = settings.show_icon;
  } else {
    next_dst_change = undefined;
  }
}
```


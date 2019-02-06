# Created using julia v1.0.2
using Dates

function genUUID()
	return parse(Int, string(Dates.value(Dates.Time(Dates.now())))[1:6])
end

function genMeIcalGoogle(fromDate, fromNumber, tillNumber)
    ret="BEGIN:VCALENDAR\nX-WR-TIMEZONE:Europe/Budapest\n"
    het1 = Week(1)
    for i in fromNumber:tillNumber
        nh = fromDate + (i-1)*het1
        np = nh + Day(5)
        ret *= "BEGIN:VEVENT\nDTSTART;VALUE=DATE:" * Dates.format(nh,"yyyymmdd") * "\n"
        ret *= "DTEND;VALUE=DATE:" * Dates.format(np,"yyyymmdd") * "\n"
        ret *= "DTSTAMP:20180904T191528Z\n"
        ret *= "TRANSP:TRANSPARENT\n"
        ret *= "SEQUENCE:0\n"
        ret *= "UID:xyz$(genUUID())$i@example.com\n"
        ret *= "DESCRIPTION:\n"
        ret *= "SUMMARY:" * "$i.hét\n"
        ret *= "PRIORITY:5\nCLASS:PUBLIC\n"
        ret *= "END:VEVENT\n"
    end
	ret *= "END:VCALENDAR"
    return ret
end

function writeToFile(name, content)
	open(name,"w") do f
	    write(f,content)
	end
end

hetfo = Date(2019,2,4)
hetfo2 = Date(2019,2,11)
vege = genMeIcalGoogle(hetfo2, 7, 14);

writeToFile("example.ics", vege)

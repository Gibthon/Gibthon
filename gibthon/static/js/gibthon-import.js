var message = function(channel, data) {
	var ifr = document.createElement('iframe');
	document.body.appendChild(ifr);
	ifr.style.display = "none";
	var frm = document.createElement('form');
	var target = 'GCD';
	ifr.contentWindow.name = target;
	frm.target = target;
	frm.action = "http://async-message-passer.appspot.com/submit?channel_name="+channel;
	frm.method = "post";
	var idt = document.createElement('textarea');
	idt.name = "content";
	idt.value = data;
	var icl = document.createElement('input');
	icl.name = "clear_channel";
	icl.value = 1;
	frm.appendChild(idt);
	frm.appendChild(icl);
	ifr.appendChild(frm);
	frm.submit();
	
	//xhr = new XMLHttpRequest();
	//xhr.open("POST", "http://async-message-passer.appspot.com/submit?channel_name="+channel, true);
	//xhr.send(data);
}
if (/partsregistry/.test(window.location.host)) {
	if (/Part\:/.test(window.location.pathname)) {
		m = /Part\:([\w_]+)/.exec(window.location.pathname);
		if (m.length == 2) {
			part = m[1];
			channel = prompt("Adding part " + part + " to Gibthon\n\nPlease enter your channel name:", "GCD-<your gibthon username here>");
			if (channel!=null && channel!=""){
				message(channel, '{"source":"pr","type":"fr","data":{"id":"'+part+'"}}');
			}
		} else {
			alert("Not found");
		}
	} else {
		alert("Not on a valid parts page");
	}
} else {
	alert("Not on a valid site");
}
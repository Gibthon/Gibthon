#!/bin/bash

#flush existing rules
iptables -F

#allow existing connections
iptables -A INPUT -m conntrack --ctstate ESTABLISHED,RELATED -j ACCEPT

#allow ssh
iptables -A INPUT -p tcp --dport ssh -j ACCEPT

#allow http
iptables -A INPUT -p tcp --dport 80 -j ACCEPT

#allow outgoing connections
iptables -P OUTPUT ACCEPT

#reject everything else
iptables -A INPUT -j REJECT


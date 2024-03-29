#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

#define N 1000000

int i,j,k,n,m,t,a[N+50],l[N+50];

map<int,int> lst;

int main(){
	ios::sync_with_stdio(0); cin.tie(0);
	cin>>n>>t;
	l[i]=0;
	for(i=1;i<=n;i++){
		cin>>a[i]; l[i]=l[i-1];
		if(lst.count(a[i])){
			if(lst[a[i]]<i-1)l[i]=max(l[i],lst[a[i]]);
		}
		if(i>=3&&a[i]==a[i-2])l[i]=max(l[i],i-2);
		lst[a[i]]=i;
	}
	while(t--){
		cin>>i>>j;
		if(j>2&&l[j]>=i)cout<<"YES\n";
		else cout<<"NO\n";
	}
}

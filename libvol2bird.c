


<!DOCTYPE html>
<html lang="en" class="">
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# object: http://ogp.me/ns/object# article: http://ogp.me/ns/article# profile: http://ogp.me/ns/profile#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta http-equiv="Content-Language" content="en">
    
    
    <title>enram/libvol2bird.c at master · NLeSC/enram</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub">
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub">
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png">
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png">
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png">
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png">
    <meta property="fb:app_id" content="1401488693436528">

      <meta content="@github" name="twitter:site" /><meta content="summary" name="twitter:card" /><meta content="NLeSC/enram" name="twitter:title" /><meta content="Contribute to enram development by creating an account on GitHub." name="twitter:description" /><meta content="https://avatars0.githubusercontent.com/u/1705862?v=3&amp;s=400" name="twitter:image:src" />
      <meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="https://avatars0.githubusercontent.com/u/1705862?v=3&amp;s=400" property="og:image" /><meta content="NLeSC/enram" property="og:title" /><meta content="https://github.com/NLeSC/enram" property="og:url" /><meta content="Contribute to enram development by creating an account on GitHub." property="og:description" />
      <meta name="browser-stats-url" content="/_stats">
    <link rel="assets" href="https://assets-cdn.github.com/">
    <link rel="xhr-socket" href="/_sockets">
    <meta name="pjax-timeout" content="1000">
    <link rel="sudo-modal" href="/sessions/sudo_modal">

    <meta name="msapplication-TileImage" content="/windows-tile.png">
    <meta name="msapplication-TileColor" content="#ffffff">
    <meta name="selected-link" value="repo_source" data-pjax-transient>
      <meta name="google-analytics" content="UA-3769691-2">

    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="9112A843:324C:1F935E8:55129521" name="octolytics-dimension-request_id" /><meta content="8268415" name="octolytics-actor-id" /><meta content="adokter" name="octolytics-actor-login" /><meta content="2c14627cf939aef86312b3ae72bcbf4993ed5a5f8d60552136254f03082d0c67" name="octolytics-actor-hash" />
    
    <meta content="Rails, view, blob#show" name="analytics-event" />

    
    <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico">


    <meta content="authenticity_token" name="csrf-param" />
<meta content="3Od48ZvKKUyEuY55ErhAxiGz9CG7jSpC2lOrVhlK0xSzItzBjljs2FOiUFhRrpE1+sNugjrX9MReDlhEGiA+1w==" name="csrf-token" />

    <link href="https://assets-cdn.github.com/assets/github-5173106e1ce2269aaefb584ff7105d11a3526c861af3895f22d7f177ed9660c2.css" media="all" rel="stylesheet" />
    <link href="https://assets-cdn.github.com/assets/github2-40d9bf14363443ccf64a2b71208f59e8739d6288d962feba121227e0608772f3.css" media="all" rel="stylesheet" />
    
    


    <meta http-equiv="x-pjax-version" content="b3746723eedf3fc96ecd834e81bbcc04">

      
  <meta name="description" content="Contribute to enram development by creating an account on GitHub.">
  <meta name="go-import" content="github.com/NLeSC/enram git https://github.com/NLeSC/enram.git">

  <meta content="1705862" name="octolytics-dimension-user_id" /><meta content="NLeSC" name="octolytics-dimension-user_login" /><meta content="14423296" name="octolytics-dimension-repository_id" /><meta content="NLeSC/enram" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="14423296" name="octolytics-dimension-repository_network_root_id" /><meta content="NLeSC/enram" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/NLeSC/enram/commits/master.atom" rel="alternate" title="Recent Commits to enram:master" type="application/atom+xml">

  </head>


  <body class="logged_in  env-production macintosh vis-public page-blob">
    <a href="#start-of-content" tabindex="1" class="accessibility-aid js-skip-to-content">Skip to content</a>
    <div class="wrapper">
      
      
      


        <div class="header header-logged-in true" role="banner">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/" data-hotkey="g d" aria-label="Homepage" data-ga-click="Header, go to dashboard, icon:logo">
  <span class="mega-octicon octicon-mark-github"></span>
</a>


      <div class="site-search repo-scope js-site-search" role="search">
          <form accept-charset="UTF-8" action="/NLeSC/enram/search" class="js-site-search-form" data-global-search-url="/search" data-repo-search-url="/NLeSC/enram/search" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
  <input type="text"
    class="js-site-search-field is-clearable"
    data-hotkey="s"
    name="q"
    placeholder="Search"
    data-global-scope-placeholder="Search GitHub"
    data-repo-scope-placeholder="Search"
    tabindex="1"
    autocapitalize="off">
  <div class="scope-badge">This repository</div>
</form>
      </div>
      <ul class="header-nav left" role="navigation">
        <li class="header-nav-item explore">
          <a class="header-nav-link" href="/explore" data-ga-click="Header, go to explore, text:explore">Explore</a>
        </li>
          <li class="header-nav-item">
            <a class="header-nav-link" href="https://gist.github.com" data-ga-click="Header, go to gist, text:gist">Gist</a>
          </li>
          <li class="header-nav-item">
            <a class="header-nav-link" href="/blog" data-ga-click="Header, go to blog, text:blog">Blog</a>
          </li>
        <li class="header-nav-item">
          <a class="header-nav-link" href="https://help.github.com" data-ga-click="Header, go to help, text:help">Help</a>
        </li>
      </ul>

    
<ul class="header-nav user-nav right" id="user-links">
  <li class="header-nav-item dropdown js-menu-container">
    <a class="header-nav-link name" href="/adokter" data-ga-click="Header, go to profile, text:username">
      <img alt="@adokter" class="avatar" data-user="8268415" height="20" src="https://avatars1.githubusercontent.com/u/8268415?v=3&amp;s=40" width="20" />
      <span class="css-truncate">
        <span class="css-truncate-target">adokter</span>
      </span>
    </a>
  </li>

  <li class="header-nav-item dropdown js-menu-container">
    <a class="header-nav-link js-menu-target tooltipped tooltipped-s" href="#" aria-label="Create new..." data-ga-click="Header, create new, icon:add">
      <span class="octicon octicon-plus"></span>
      <span class="dropdown-caret"></span>
    </a>

    <div class="dropdown-menu-content js-menu-content">
      
<ul class="dropdown-menu">
  <li>
    <a href="/new" data-ga-click="Header, create new repository, icon:repo"><span class="octicon octicon-repo"></span> New repository</a>
  </li>
  <li>
    <a href="/organizations/new" data-ga-click="Header, create new organization, icon:organization"><span class="octicon octicon-organization"></span> New organization</a>
  </li>


    <li class="dropdown-divider"></li>
    <li class="dropdown-header">
      <span title="NLeSC/enram">This repository</span>
    </li>
      <li>
        <a href="/NLeSC/enram/issues/new" data-ga-click="Header, create new issue, icon:issue"><span class="octicon octicon-issue-opened"></span> New issue</a>
      </li>
</ul>

    </div>
  </li>

  <li class="header-nav-item">
        <a href="/NLeSC/enram/notifications" aria-label="You have unread notifications in this repository" class="header-nav-link notification-indicator tooltipped tooltipped-s" data-ga-click="Header, go to notifications, icon:unread" data-hotkey="g n">
        <span class="mail-status unread"></span>
        <span class="octicon octicon-inbox"></span>
</a>
  </li>

  <li class="header-nav-item">
    <a class="header-nav-link tooltipped tooltipped-s" href="/settings/profile" id="account_settings" aria-label="Settings" data-ga-click="Header, go to settings, icon:settings">
      <span class="octicon octicon-gear"></span>
    </a>
  </li>

  <li class="header-nav-item">
    <form accept-charset="UTF-8" action="/logout" class="logout-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="9oKrnLg9FqWobFAhnASXVnt7viX15hkTqpoLDOt5ns1bqxUGp6epFCl5wgtBkCAhkNukk4wozO3rbv5Y4OLqVA==" /></div>
      <button class="header-nav-link sign-out-button tooltipped tooltipped-s" aria-label="Sign out" data-ga-click="Header, sign out, icon:logout">
        <span class="octicon octicon-sign-out"></span>
      </button>
</form>  </li>

</ul>


    
  </div>
</div>

        

        


      <div id="start-of-content" class="accessibility-aid"></div>
          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    <div id="js-flash-container">
      
    </div>
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        
<ul class="pagehead-actions">

  <li>
      <form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="zB8dipnLQNu3MHq9CWRFSUmCCN0eJPjgaIZCHDgUlXKJgpYkl2Pr8P6ngxThYa8r34v6YCGdVEbOKLmDaaTRFA==" /></div>    <input id="repository_id" name="repository_id" type="hidden" value="14423296" />

      <div class="select-menu js-menu-container js-select-menu">
        <a href="/NLeSC/enram/subscription"
          class="btn btn-sm btn-with-count select-menu-button js-menu-target" role="button" tabindex="0" aria-haspopup="true"
          data-ga-click="Repository, click Watch settings, action:blob#show">
          <span class="js-select-button">
            <span class="octicon octicon-eye"></span>
            Watch
          </span>
        </a>
        <a class="social-count js-social-count" href="/NLeSC/enram/watchers">
          0
        </a>

        <div class="select-menu-modal-holder">
          <div class="select-menu-modal subscription-menu-modal js-menu-content" aria-hidden="true">
            <div class="select-menu-header">
              <span class="select-menu-title">Notifications</span>
              <span class="octicon octicon-x js-menu-close" role="button" aria-label="Close"></span>
            </div>

            <div class="select-menu-list js-navigation-container" role="menu">

              <div class="select-menu-item js-navigation-item selected" role="menuitem" tabindex="0">
                <span class="select-menu-item-icon octicon octicon-check"></span>
                <div class="select-menu-item-text">
                  <input checked="checked" id="do_included" name="do" type="radio" value="included" />
                  <span class="select-menu-item-heading">Not watching</span>
                  <span class="description">Be notified when participating or @mentioned.</span>
                  <span class="js-select-button-text hidden-select-button-text">
                    <span class="octicon octicon-eye"></span>
                    Watch
                  </span>
                </div>
              </div>

              <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
                <span class="select-menu-item-icon octicon octicon octicon-check"></span>
                <div class="select-menu-item-text">
                  <input id="do_subscribed" name="do" type="radio" value="subscribed" />
                  <span class="select-menu-item-heading">Watching</span>
                  <span class="description">Be notified of all conversations.</span>
                  <span class="js-select-button-text hidden-select-button-text">
                    <span class="octicon octicon-eye"></span>
                    Unwatch
                  </span>
                </div>
              </div>

              <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
                <span class="select-menu-item-icon octicon octicon-check"></span>
                <div class="select-menu-item-text">
                  <input id="do_ignore" name="do" type="radio" value="ignore" />
                  <span class="select-menu-item-heading">Ignoring</span>
                  <span class="description">Never be notified.</span>
                  <span class="js-select-button-text hidden-select-button-text">
                    <span class="octicon octicon-mute"></span>
                    Stop ignoring
                  </span>
                </div>
              </div>

            </div>

          </div>
        </div>
      </div>
</form>
  </li>

  <li>
    
  <div class="js-toggler-container js-social-container starring-container ">

    <form accept-charset="UTF-8" action="/NLeSC/enram/unstar" class="js-toggler-form starred js-unstar-button" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="VZBWojCG4W/UPr79zsjKDpPWS6wsFh0K2F3LWbxJz8OExN+u9ltJoIYSK0DC2RXPhDqKpU3OfM+XxkI2EffyMg==" /></div>
      <button
        class="btn btn-sm btn-with-count js-toggler-target"
        aria-label="Unstar this repository" title="Unstar NLeSC/enram"
        data-ga-click="Repository, click unstar button, action:blob#show; text:Unstar">
        <span class="octicon octicon-star"></span>
        Unstar
      </button>
        <a class="social-count js-social-count" href="/NLeSC/enram/stargazers">
          0
        </a>
</form>
    <form accept-charset="UTF-8" action="/NLeSC/enram/star" class="js-toggler-form unstarred js-star-button" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="inwdey7xPiIThNNBTBRK9tR0RyPNi+fMoWYLo+PHPNzsozf7SLvVTSh6q8dV8foaNfE/7Dr/his76X0otJ0+qQ==" /></div>
      <button
        class="btn btn-sm btn-with-count js-toggler-target"
        aria-label="Star this repository" title="Star NLeSC/enram"
        data-ga-click="Repository, click star button, action:blob#show; text:Star">
        <span class="octicon octicon-star"></span>
        Star
      </button>
        <a class="social-count js-social-count" href="/NLeSC/enram/stargazers">
          0
        </a>
</form>  </div>

  </li>

        <li>
          <form accept-charset="UTF-8" action="/NLeSC/enram/fork" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="RK9w2/JndLxc6MW0Iv3ZYAB5hvFueIH6MRiR+d5Db+pbD+cUYgb0mIP2EnenC6TyhxljTNM5KxDndyJkbx1ghg==" /></div>
            <button
                type="submit"
                class="btn btn-sm btn-with-count"
                data-ga-click="Repository, show fork modal, action:blob#show; text:Fork"
                title="Fork your own copy of NLeSC/enram to your account"
                aria-label="Fork your own copy of NLeSC/enram to your account">
              <span class="octicon octicon-repo-forked"></span>
              Fork
            </button>
            <a href="/NLeSC/enram/network" class="social-count">1</a>
</form>        </li>

</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="mega-octicon octicon-repo"></span>
          <span class="author"><a href="/NLeSC" class="url fn" itemprop="url" rel="author"><span itemprop="title">NLeSC</span></a></span><!--
       --><span class="path-divider">/</span><!--
       --><strong><a href="/NLeSC/enram" class="js-current-repository" data-pjax="#js-repo-pjax-container">enram</a></strong>

          <span class="page-context-loader">
            <img alt="" height="16" src="https://assets-cdn.github.com/assets/spinners/octocat-spinner-32-e513294efa576953719e4e2de888dd9cf929b7d62ed8d05f25e731d02452ab6c.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">
      <div class="repository-with-sidebar repo-container new-discussion-timeline  ">
        <div class="repository-sidebar clearfix">
            
<nav class="sunken-menu repo-nav js-repo-nav js-sidenav-container-pjax js-octicon-loaders"
     role="navigation"
     data-pjax="#js-repo-pjax-container"
     data-issue-count-url="/NLeSC/enram/issues/counts">
  <ul class="sunken-menu-group">
    <li class="tooltipped tooltipped-w" aria-label="Code">
      <a href="/NLeSC/enram" aria-label="Code" class="selected js-selected-navigation-item sunken-menu-item" data-hotkey="g c" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /NLeSC/enram">
        <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
        <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/assets/spinners/octocat-spinner-32-e513294efa576953719e4e2de888dd9cf929b7d62ed8d05f25e731d02452ab6c.gif" width="16" />
</a>    </li>

      <li class="tooltipped tooltipped-w" aria-label="Issues">
        <a href="/NLeSC/enram/issues" aria-label="Issues" class="js-selected-navigation-item sunken-menu-item" data-hotkey="g i" data-selected-links="repo_issues repo_labels repo_milestones /NLeSC/enram/issues">
          <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
          <span class="js-issue-replace-counter"></span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/assets/spinners/octocat-spinner-32-e513294efa576953719e4e2de888dd9cf929b7d62ed8d05f25e731d02452ab6c.gif" width="16" />
</a>      </li>

    <li class="tooltipped tooltipped-w" aria-label="Pull requests">
      <a href="/NLeSC/enram/pulls" aria-label="Pull requests" class="js-selected-navigation-item sunken-menu-item" data-hotkey="g p" data-selected-links="repo_pulls /NLeSC/enram/pulls">
          <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull requests</span>
          <span class="js-pull-replace-counter"></span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/assets/spinners/octocat-spinner-32-e513294efa576953719e4e2de888dd9cf929b7d62ed8d05f25e731d02452ab6c.gif" width="16" />
</a>    </li>


      <li class="tooltipped tooltipped-w" aria-label="Wiki">
        <a href="/NLeSC/enram/wiki" aria-label="Wiki" class="js-selected-navigation-item sunken-menu-item" data-hotkey="g w" data-selected-links="repo_wiki /NLeSC/enram/wiki">
          <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/assets/spinners/octocat-spinner-32-e513294efa576953719e4e2de888dd9cf929b7d62ed8d05f25e731d02452ab6c.gif" width="16" />
</a>      </li>
  </ul>
  <div class="sunken-menu-separator"></div>
  <ul class="sunken-menu-group">

    <li class="tooltipped tooltipped-w" aria-label="Pulse">
      <a href="/NLeSC/enram/pulse" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-selected-links="pulse /NLeSC/enram/pulse">
        <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
        <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/assets/spinners/octocat-spinner-32-e513294efa576953719e4e2de888dd9cf929b7d62ed8d05f25e731d02452ab6c.gif" width="16" />
</a>    </li>

    <li class="tooltipped tooltipped-w" aria-label="Graphs">
      <a href="/NLeSC/enram/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-selected-links="repo_graphs repo_contributors /NLeSC/enram/graphs">
        <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
        <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/assets/spinners/octocat-spinner-32-e513294efa576953719e4e2de888dd9cf929b7d62ed8d05f25e731d02452ab6c.gif" width="16" />
</a>    </li>
  </ul>


</nav>

              <div class="only-with-full-nav">
                  
<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=clone">
  <h3><span class="text-emphasized">HTTPS</span> clone URL</h3>
  <div class="input-group js-zeroclipboard-container">
    <input type="text" class="input-mini input-monospace js-url-field js-zeroclipboard-target"
           value="https://github.com/NLeSC/enram.git" readonly="readonly">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard btn btn-sm zeroclipboard-button" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  
<div class="clone-url "
  data-protocol-type="ssh"
  data-url="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=clone">
  <h3><span class="text-emphasized">SSH</span> clone URL</h3>
  <div class="input-group js-zeroclipboard-container">
    <input type="text" class="input-mini input-monospace js-url-field js-zeroclipboard-target"
           value="git@github.com:NLeSC/enram.git" readonly="readonly">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard btn btn-sm zeroclipboard-button" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  
<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=clone">
  <h3><span class="text-emphasized">Subversion</span> checkout URL</h3>
  <div class="input-group js-zeroclipboard-container">
    <input type="text" class="input-mini input-monospace js-url-field js-zeroclipboard-target"
           value="https://github.com/NLeSC/enram" readonly="readonly">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard btn btn-sm zeroclipboard-button" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>



<p class="clone-options">You can clone with
  <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>, <a href="#" class="js-clone-selector" data-protocol="ssh">SSH</a>, or <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>.
  <a href="https://help.github.com/articles/which-remote-url-should-i-use" class="help tooltipped tooltipped-n" aria-label="Get help on which URL is right for you.">
    <span class="octicon octicon-question"></span>
  </a>
</p>

  <a href="github-mac://openRepo/https://github.com/NLeSC/enram" class="btn btn-sm sidebar-button" title="Save NLeSC/enram to your computer and use it in GitHub Desktop." aria-label="Save NLeSC/enram to your computer and use it in GitHub Desktop.">
    <span class="octicon octicon-device-desktop"></span>
    Clone in Desktop
  </a>


                <a href="/NLeSC/enram/archive/master.zip"
                   class="btn btn-sm sidebar-button"
                   aria-label="Download the contents of NLeSC/enram as a zip file"
                   title="Download the contents of NLeSC/enram as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
              </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          

<a href="/NLeSC/enram/blob/b93e42d608e44117331f5e2743d73962234d8a39/ncradar/lib/vol2bird/libvol2bird.c" class="hidden js-permalink-shortcut" data-hotkey="y">Permalink</a>

<!-- blob contrib key: blob_contributors:v21:4bd73b9121ab47e790e3251e88cbae1c -->

<div class="file-navigation js-zeroclipboard-container">
  
<div class="select-menu js-menu-container js-select-menu left">
  <span class="btn btn-sm select-menu-button js-menu-target css-truncate" data-hotkey="w"
    data-master-branch="master"
    data-ref="master"
    title="master"
    role="button" aria-label="Switch branches or tags" tabindex="0" aria-haspopup="true">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button css-truncate-target">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax aria-hidden="true">

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-x js-menu-close" role="button" aria-label="Close"></span>
      </div>

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Filter branches/tags" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" data-filter-placeholder="Filter branches/tags" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" data-filter-placeholder="Find a tag…" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/NLeSC/enram/blob/develop/ncradar/lib/vol2bird/libvol2bird.c"
               data-name="develop"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="develop">
                develop
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/NLeSC/enram/blob/issue70/ncradar/lib/vol2bird/libvol2bird.c"
               data-name="issue70"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="issue70">
                issue70
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/NLeSC/enram/blob/issue89/ncradar/lib/vol2bird/libvol2bird.c"
               data-name="issue89"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="issue89">
                issue89
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/NLeSC/enram/blob/javascript-visualization/ncradar/lib/vol2bird/libvol2bird.c"
               data-name="javascript-visualization"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="javascript-visualization">
                javascript-visualization
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open selected"
               href="/NLeSC/enram/blob/master/ncradar/lib/vol2bird/libvol2bird.c"
               data-name="master"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="master">
                master
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/NLeSC/enram/blob/profileDataGetter/ncradar/lib/vol2bird/libvol2bird.c"
               data-name="profileDataGetter"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="profileDataGetter">
                profileDataGetter
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/NLeSC/enram/blob/withOrWithoutSvdfit/ncradar/lib/vol2bird/libvol2bird.c"
               data-name="withOrWithoutSvdfit"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="withOrWithoutSvdfit">
                withOrWithoutSvdfit
              </span>
            </a>
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div>

    </div>
  </div>
</div>

  <div class="btn-group right">
    <a href="/NLeSC/enram/find/master"
          class="js-show-file-finder btn btn-sm empty-icon tooltipped tooltipped-s"
          data-pjax
          data-hotkey="t"
          aria-label="Quickly jump between files">
      <span class="octicon octicon-list-unordered"></span>
    </a>
    <button aria-label="Copy file path to clipboard" class="js-zeroclipboard btn btn-sm zeroclipboard-button" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
  </div>

  <div class="breadcrumb js-zeroclipboard-target">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/NLeSC/enram" class="" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">enram</span></a></span></span><span class="separator">/</span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/NLeSC/enram/tree/master/ncradar" class="" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">ncradar</span></a></span><span class="separator">/</span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/NLeSC/enram/tree/master/ncradar/lib" class="" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">lib</span></a></span><span class="separator">/</span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/NLeSC/enram/tree/master/ncradar/lib/vol2bird" class="" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">vol2bird</span></a></span><span class="separator">/</span><strong class="final-path">libvol2bird.c</strong>
  </div>
</div>


  <div class="commit file-history-tease">
    <div class="file-history-tease-header">
        <img alt="@jspaaks" class="avatar" data-user="4558105" height="24" src="https://avatars2.githubusercontent.com/u/4558105?v=3&amp;s=48" width="24" />
        <span class="author"><a href="/jspaaks" rel="contributor">jspaaks</a></span>
        <time datetime="2015-03-09T13:54:53Z" is="relative-time">Mar 9, 2015</time>
        <div class="commit-title">
            <a href="/NLeSC/enram/commit/b93e42d608e44117331f5e2743d73962234d8a39" class="message" data-pjax="true" title="set a couple of printer methods to static and moved them to the .c file where they belong">set a couple of printer methods to static and moved them to the .c fi…</a>
        </div>
    </div>

    <div class="participation">
      <p class="quickstat">
        <a href="#blob_contributors_box" rel="facebox">
          <strong>1</strong>
           contributor
        </a>
      </p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list">
          <li class="facebox-user-list-item">
            <img alt="@jspaaks" data-user="4558105" height="24" src="https://avatars2.githubusercontent.com/u/4558105?v=3&amp;s=48" width="24" />
            <a href="/jspaaks">jspaaks</a>
          </li>
      </ul>
    </div>
  </div>

<div class="file">
  <div class="file-header">
    <div class="file-actions">

      <div class="btn-group">
        <a href="/NLeSC/enram/raw/master/ncradar/lib/vol2bird/libvol2bird.c" class="btn btn-sm " id="raw-url">Raw</a>
          <a href="/NLeSC/enram/blame/master/ncradar/lib/vol2bird/libvol2bird.c" class="btn btn-sm js-update-url-with-hash">Blame</a>
        <a href="/NLeSC/enram/commits/master/ncradar/lib/vol2bird/libvol2bird.c" class="btn btn-sm " rel="nofollow">History</a>
      </div>

        <a class="octicon-btn tooltipped tooltipped-nw"
           href="github-mac://openRepo/https://github.com/NLeSC/enram?branch=master&amp;filepath=ncradar%2Flib%2Fvol2bird%2Flibvol2bird.c"
           aria-label="Open this file in GitHub for Mac">
            <span class="octicon octicon-device-desktop"></span>
        </a>

            <form accept-charset="UTF-8" action="/NLeSC/enram/edit/master/ncradar/lib/vol2bird/libvol2bird.c" class="inline-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="ZljpYSlwG/py524Wci9dJ656OM6g6i3/XU89fG6tGNmyyrqnUCJxN9EqXZA6TVwkyQ4aHSlLWeS4gARi/HySJw==" /></div>
              <button class="octicon-btn tooltipped tooltipped-n" type="submit" aria-label="Clicking this button will fork this project so you can edit the file" data-hotkey="e" data-disable-with>
                <span class="octicon octicon-pencil"></span>
              </button>
</form>
          <form accept-charset="UTF-8" action="/NLeSC/enram/delete/master/ncradar/lib/vol2bird/libvol2bird.c" class="inline-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="gkSm6YPyBKpZVUau8NAA2s/sxKMEcwF8SGGHIKamW5KfwjTJ0TE6EFmNBPYkNVIOiRQZwofDUhr5e8/keNrr3g==" /></div>
            <button class="octicon-btn octicon-btn-danger tooltipped tooltipped-n" type="submit" aria-label="Fork this project and delete file" data-disable-with>
              <span class="octicon octicon-trashcan"></span>
            </button>
</form>    </div>

    <div class="file-info">
        3133 lines (2262 sloc)
        <span class="file-info-divider"></span>
      104.774 kb
    </div>
  </div>
  
  <div class="blob-wrapper data type-c">
      <table class="highlight tab-size-8 js-file-line-container">
      <tr>
        <td id="L1" class="blob-num js-line-number" data-line-number="1"></td>
        <td id="LC1" class="blob-code js-file-line"><span class="pl-c">/*</span></td>
      </tr>
      <tr>
        <td id="L2" class="blob-num js-line-number" data-line-number="2"></td>
        <td id="LC2" class="blob-code js-file-line"><span class="pl-c"> * Copyright 2013 Netherlands eScience Center</span></td>
      </tr>
      <tr>
        <td id="L3" class="blob-num js-line-number" data-line-number="3"></td>
        <td id="LC3" class="blob-code js-file-line"><span class="pl-c"> *</span></td>
      </tr>
      <tr>
        <td id="L4" class="blob-num js-line-number" data-line-number="4"></td>
        <td id="LC4" class="blob-code js-file-line"><span class="pl-c"> * Licensed under the Apache License, Version 2.0 (the &quot;License&quot;);</span></td>
      </tr>
      <tr>
        <td id="L5" class="blob-num js-line-number" data-line-number="5"></td>
        <td id="LC5" class="blob-code js-file-line"><span class="pl-c"> * you may not use this file except in compliance with the License.</span></td>
      </tr>
      <tr>
        <td id="L6" class="blob-num js-line-number" data-line-number="6"></td>
        <td id="LC6" class="blob-code js-file-line"><span class="pl-c"> * You may obtain a copy of the License at</span></td>
      </tr>
      <tr>
        <td id="L7" class="blob-num js-line-number" data-line-number="7"></td>
        <td id="LC7" class="blob-code js-file-line"><span class="pl-c"> *</span></td>
      </tr>
      <tr>
        <td id="L8" class="blob-num js-line-number" data-line-number="8"></td>
        <td id="LC8" class="blob-code js-file-line"><span class="pl-c"> * http://www.apache.org/licenses/LICENSE-2.0</span></td>
      </tr>
      <tr>
        <td id="L9" class="blob-num js-line-number" data-line-number="9"></td>
        <td id="LC9" class="blob-code js-file-line"><span class="pl-c"> *</span></td>
      </tr>
      <tr>
        <td id="L10" class="blob-num js-line-number" data-line-number="10"></td>
        <td id="LC10" class="blob-code js-file-line"><span class="pl-c"> * Unless required by applicable law or agreed to in writing, software</span></td>
      </tr>
      <tr>
        <td id="L11" class="blob-num js-line-number" data-line-number="11"></td>
        <td id="LC11" class="blob-code js-file-line"><span class="pl-c"> * distributed under the License is distributed on an &quot;AS IS&quot; BASIS,</span></td>
      </tr>
      <tr>
        <td id="L12" class="blob-num js-line-number" data-line-number="12"></td>
        <td id="LC12" class="blob-code js-file-line"><span class="pl-c"> * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.</span></td>
      </tr>
      <tr>
        <td id="L13" class="blob-num js-line-number" data-line-number="13"></td>
        <td id="LC13" class="blob-code js-file-line"><span class="pl-c"> * See the License for the specific language governing permissions and</span></td>
      </tr>
      <tr>
        <td id="L14" class="blob-num js-line-number" data-line-number="14"></td>
        <td id="LC14" class="blob-code js-file-line"><span class="pl-c"> * limitations under the License.</span></td>
      </tr>
      <tr>
        <td id="L15" class="blob-num js-line-number" data-line-number="15"></td>
        <td id="LC15" class="blob-code js-file-line"><span class="pl-c"> */</span></td>
      </tr>
      <tr>
        <td id="L16" class="blob-num js-line-number" data-line-number="16"></td>
        <td id="LC16" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L17" class="blob-num js-line-number" data-line-number="17"></td>
        <td id="LC17" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L18" class="blob-num js-line-number" data-line-number="18"></td>
        <td id="LC18" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>stdio.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L19" class="blob-num js-line-number" data-line-number="19"></td>
        <td id="LC19" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>confuse.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L20" class="blob-num js-line-number" data-line-number="20"></td>
        <td id="LC20" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>stdlib.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L21" class="blob-num js-line-number" data-line-number="21"></td>
        <td id="LC21" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>math.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L22" class="blob-num js-line-number" data-line-number="22"></td>
        <td id="LC22" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>vertical_profile.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L23" class="blob-num js-line-number" data-line-number="23"></td>
        <td id="LC23" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L24" class="blob-num js-line-number" data-line-number="24"></td>
        <td id="LC24" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>polarvolume.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L25" class="blob-num js-line-number" data-line-number="25"></td>
        <td id="LC25" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>polarscan.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L26" class="blob-num js-line-number" data-line-number="26"></td>
        <td id="LC26" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>libvol2bird.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L27" class="blob-num js-line-number" data-line-number="27"></td>
        <td id="LC27" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>libsvdfit.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L28" class="blob-num js-line-number" data-line-number="28"></td>
        <td id="LC28" class="blob-code js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>constants.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L29" class="blob-num js-line-number" data-line-number="29"></td>
        <td id="LC29" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L30" class="blob-num js-line-number" data-line-number="30"></td>
        <td id="LC30" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L31" class="blob-num js-line-number" data-line-number="31"></td>
        <td id="LC31" class="blob-code js-file-line"><span class="pl-c">// non-public function prototypes (local to this file/translation unit)</span></td>
      </tr>
      <tr>
        <td id="L32" class="blob-num js-line-number" data-line-number="32"></td>
        <td id="LC32" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L33" class="blob-num js-line-number" data-line-number="33"></td>
        <td id="LC33" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">analyzeCells</span>(<span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *dbzImage, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *vradImage,</td>
      </tr>
      <tr>
        <td id="L34" class="blob-num js-line-number" data-line-number="34"></td>
        <td id="LC34" class="blob-code js-file-line">                        <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *texImage, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *clutterImage, <span class="pl-k">int</span> *cellImage,</td>
      </tr>
      <tr>
        <td id="L35" class="blob-num js-line-number" data-line-number="35"></td>
        <td id="LC35" class="blob-code js-file-line">                        <span class="pl-k">const</span> SCANMETA *dbzMeta, <span class="pl-k">const</span> SCANMETA *vradMeta, <span class="pl-k">const</span> SCANMETA *texMeta, <span class="pl-k">const</span> SCANMETA *clutterMeta,</td>
      </tr>
      <tr>
        <td id="L36" class="blob-num js-line-number" data-line-number="36"></td>
        <td id="LC36" class="blob-code js-file-line">                        <span class="pl-k">const</span> <span class="pl-k">int</span> nCells, <span class="pl-k">const</span> <span class="pl-k">int</span> useStaticClutterData);</td>
      </tr>
      <tr>
        <td id="L37" class="blob-num js-line-number" data-line-number="37"></td>
        <td id="LC37" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L38" class="blob-num js-line-number" data-line-number="38"></td>
        <td id="LC38" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> <span class="pl-en">calcDist</span>(<span class="pl-k">const</span> <span class="pl-k">int</span> range1, <span class="pl-k">const</span> <span class="pl-k">int</span> azim1, <span class="pl-k">const</span> <span class="pl-k">int</span> range2, <span class="pl-k">const</span> <span class="pl-k">int</span> azim2, <span class="pl-k">const</span> <span class="pl-k">float</span> rscale, <span class="pl-k">const</span> <span class="pl-k">float</span> ascale);</td>
      </tr>
      <tr>
        <td id="L39" class="blob-num js-line-number" data-line-number="39"></td>
        <td id="LC39" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L40" class="blob-num js-line-number" data-line-number="40"></td>
        <td id="LC40" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">calcTexture</span>(<span class="pl-k">unsigned</span> <span class="pl-k">char</span> *texImage, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *vradImage, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *dbzImage,</td>
      </tr>
      <tr>
        <td id="L41" class="blob-num js-line-number" data-line-number="41"></td>
        <td id="LC41" class="blob-code js-file-line">                        <span class="pl-k">const</span> SCANMETA *texMeta, <span class="pl-k">const</span> SCANMETA *vradMeta, <span class="pl-k">const</span> SCANMETA *dbzMeta);</td>
      </tr>
      <tr>
        <td id="L42" class="blob-num js-line-number" data-line-number="42"></td>
        <td id="LC42" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L43" class="blob-num js-line-number" data-line-number="43"></td>
        <td id="LC43" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">classifyGatesSimple</span>(<span class="pl-k">void</span>);</td>
      </tr>
      <tr>
        <td id="L44" class="blob-num js-line-number" data-line-number="44"></td>
        <td id="LC44" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L45" class="blob-num js-line-number" data-line-number="45"></td>
        <td id="LC45" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">constructorInt</span>(SCANMETA* meta, <span class="pl-k">int</span>* image, PolarScan_t* scan, <span class="pl-k">const</span> <span class="pl-k">int</span> nGlobal, <span class="pl-k">const</span> <span class="pl-k">int</span> initValue);</td>
      </tr>
      <tr>
        <td id="L46" class="blob-num js-line-number" data-line-number="46"></td>
        <td id="LC46" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L47" class="blob-num js-line-number" data-line-number="47"></td>
        <td id="LC47" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">constructorUChar</span>(SCANMETA* meta, <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* image, PolarScan_t* scan, <span class="pl-k">const</span> <span class="pl-k">int</span> nGlobal, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> initValue);</td>
      </tr>
      <tr>
        <td id="L48" class="blob-num js-line-number" data-line-number="48"></td>
        <td id="LC48" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L49" class="blob-num js-line-number" data-line-number="49"></td>
        <td id="LC49" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">constructPointsArray</span>(PolarVolume_t* volume);</td>
      </tr>
      <tr>
        <td id="L50" class="blob-num js-line-number" data-line-number="50"></td>
        <td id="LC50" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L51" class="blob-num js-line-number" data-line-number="51"></td>
        <td id="LC51" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">detNumberOfGates</span>(<span class="pl-k">const</span> <span class="pl-k">int</span> iLayer, <span class="pl-k">const</span> <span class="pl-k">float</span> rangeScale, <span class="pl-k">const</span> <span class="pl-k">float</span> elevAngle,</td>
      </tr>
      <tr>
        <td id="L52" class="blob-num js-line-number" data-line-number="52"></td>
        <td id="LC52" class="blob-code js-file-line">                            <span class="pl-k">const</span> <span class="pl-k">int</span> nRang, <span class="pl-k">const</span> <span class="pl-k">int</span> nAzim, <span class="pl-k">const</span> <span class="pl-k">float</span> radarHeight);</td>
      </tr>
      <tr>
        <td id="L53" class="blob-num js-line-number" data-line-number="53"></td>
        <td id="LC53" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L54" class="blob-num js-line-number" data-line-number="54"></td>
        <td id="LC54" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">detSvdfitArraySize</span>(PolarVolume_t* volume);</td>
      </tr>
      <tr>
        <td id="L55" class="blob-num js-line-number" data-line-number="55"></td>
        <td id="LC55" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L56" class="blob-num js-line-number" data-line-number="56"></td>
        <td id="LC56" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">findWeatherCells</span>(<span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *dbzImage, <span class="pl-k">int</span> *cellImage, <span class="pl-k">const</span> SCANMETA *dbzMeta);</td>
      </tr>
      <tr>
        <td id="L57" class="blob-num js-line-number" data-line-number="57"></td>
        <td id="LC57" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L58" class="blob-num js-line-number" data-line-number="58"></td>
        <td id="LC58" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">findNearbyGateIndex</span>(<span class="pl-k">const</span> <span class="pl-k">int</span> nAzimParent, <span class="pl-k">const</span> <span class="pl-k">int</span> nRangParent, <span class="pl-k">const</span> <span class="pl-k">int</span> iParent,</td>
      </tr>
      <tr>
        <td id="L59" class="blob-num js-line-number" data-line-number="59"></td>
        <td id="LC59" class="blob-code js-file-line">                               <span class="pl-k">const</span> <span class="pl-k">int</span> nAzimChild,  <span class="pl-k">const</span> <span class="pl-k">int</span> nRangChild,  <span class="pl-k">const</span> <span class="pl-k">int</span> iChild);</td>
      </tr>
      <tr>
        <td id="L60" class="blob-num js-line-number" data-line-number="60"></td>
        <td id="LC60" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L61" class="blob-num js-line-number" data-line-number="61"></td>
        <td id="LC61" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">fringeCells</span>(<span class="pl-k">int</span> *cellImage,<span class="pl-k">int</span> nRang, <span class="pl-k">int</span> nAzim, <span class="pl-k">float</span> aScale, <span class="pl-k">float</span> rScale);</td>
      </tr>
      <tr>
        <td id="L62" class="blob-num js-line-number" data-line-number="62"></td>
        <td id="LC62" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L63" class="blob-num js-line-number" data-line-number="63"></td>
        <td id="LC63" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">getListOfSelectedGates</span>(<span class="pl-k">const</span> SCANMETA* vradMeta, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *vradImage,</td>
      </tr>
      <tr>
        <td id="L64" class="blob-num js-line-number" data-line-number="64"></td>
        <td id="LC64" class="blob-code js-file-line">                                  <span class="pl-k">const</span> SCANMETA* dbzMeta, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *dbzImage,</td>
      </tr>
      <tr>
        <td id="L65" class="blob-num js-line-number" data-line-number="65"></td>
        <td id="LC65" class="blob-code js-file-line">                                  <span class="pl-k">const</span> <span class="pl-k">int</span> *cellImage,</td>
      </tr>
      <tr>
        <td id="L66" class="blob-num js-line-number" data-line-number="66"></td>
        <td id="LC66" class="blob-code js-file-line">                                  <span class="pl-k">const</span> <span class="pl-k">float</span> altitudeMin, <span class="pl-k">const</span> <span class="pl-k">float</span> altitudeMax,</td>
      </tr>
      <tr>
        <td id="L67" class="blob-num js-line-number" data-line-number="67"></td>
        <td id="LC67" class="blob-code js-file-line">                                  <span class="pl-k">float</span>* points, <span class="pl-k">int</span> iPoint);</td>
      </tr>
      <tr>
        <td id="L68" class="blob-num js-line-number" data-line-number="68"></td>
        <td id="LC68" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L69" class="blob-num js-line-number" data-line-number="69"></td>
        <td id="LC69" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">hasAzimuthGap</span>(<span class="pl-k">const</span> <span class="pl-k">float</span> *points, <span class="pl-k">const</span> <span class="pl-k">int</span> nPoints);</td>
      </tr>
      <tr>
        <td id="L70" class="blob-num js-line-number" data-line-number="70"></td>
        <td id="LC70" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L71" class="blob-num js-line-number" data-line-number="71"></td>
        <td id="LC71" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">includeGate</span>(<span class="pl-k">const</span> <span class="pl-k">int</span> iProfileType, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">int</span> gateCode);</td>
      </tr>
      <tr>
        <td id="L72" class="blob-num js-line-number" data-line-number="72"></td>
        <td id="LC72" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L73" class="blob-num js-line-number" data-line-number="73"></td>
        <td id="LC73" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">mapDataFromRave</span>(PolarScan_t* scan, SCANMETA *meta, </td>
      </tr>
      <tr>
        <td id="L74" class="blob-num js-line-number" data-line-number="74"></td>
        <td id="LC74" class="blob-code js-file-line">                           <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *values, <span class="pl-k">char</span> *paramStr);</td>
      </tr>
      <tr>
        <td id="L75" class="blob-num js-line-number" data-line-number="75"></td>
        <td id="LC75" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L76" class="blob-num js-line-number" data-line-number="76"></td>
        <td id="LC76" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">mapDataToRave</span>(<span class="pl-k">void</span>);</td>
      </tr>
      <tr>
        <td id="L77" class="blob-num js-line-number" data-line-number="77"></td>
        <td id="LC77" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L78" class="blob-num js-line-number" data-line-number="78"></td>
        <td id="LC78" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">printGateCode</span>(<span class="pl-k">char</span>* flags, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">int</span> gateCode);</td>
      </tr>
      <tr>
        <td id="L79" class="blob-num js-line-number" data-line-number="79"></td>
        <td id="LC79" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L80" class="blob-num js-line-number" data-line-number="80"></td>
        <td id="LC80" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">printImageInt</span>(<span class="pl-k">const</span> SCANMETA* meta, <span class="pl-k">const</span> <span class="pl-k">int</span>* imageInt);</td>
      </tr>
      <tr>
        <td id="L81" class="blob-num js-line-number" data-line-number="81"></td>
        <td id="LC81" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L82" class="blob-num js-line-number" data-line-number="82"></td>
        <td id="LC82" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">printImageUChar</span>(<span class="pl-k">const</span> SCANMETA* meta, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* imageUChar);</td>
      </tr>
      <tr>
        <td id="L83" class="blob-num js-line-number" data-line-number="83"></td>
        <td id="LC83" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L84" class="blob-num js-line-number" data-line-number="84"></td>
        <td id="LC84" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">printMeta</span>(<span class="pl-k">const</span> SCANMETA* meta, <span class="pl-k">const</span> <span class="pl-k">char</span>* varName);</td>
      </tr>
      <tr>
        <td id="L85" class="blob-num js-line-number" data-line-number="85"></td>
        <td id="LC85" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L86" class="blob-num js-line-number" data-line-number="86"></td>
        <td id="LC86" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">printProfile</span>(<span class="pl-k">void</span>);</td>
      </tr>
      <tr>
        <td id="L87" class="blob-num js-line-number" data-line-number="87"></td>
        <td id="LC87" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L88" class="blob-num js-line-number" data-line-number="88"></td>
        <td id="LC88" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">sortCells</span>(CELLPROP *cellProp, <span class="pl-k">const</span> <span class="pl-k">int</span> nCells);</td>
      </tr>
      <tr>
        <td id="L89" class="blob-num js-line-number" data-line-number="89"></td>
        <td id="LC89" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L90" class="blob-num js-line-number" data-line-number="90"></td>
        <td id="LC90" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">updateFlagFieldsInPointsArray</span>(<span class="pl-k">const</span> <span class="pl-k">float</span>* yObs, <span class="pl-k">const</span> <span class="pl-k">float</span>* yFitted, <span class="pl-k">const</span> <span class="pl-k">int</span>* includedIndex, </td>
      </tr>
      <tr>
        <td id="L91" class="blob-num js-line-number" data-line-number="91"></td>
        <td id="LC91" class="blob-code js-file-line">                                          <span class="pl-k">const</span> <span class="pl-k">int</span> nPointsIncluded, <span class="pl-k">float</span>* points);</td>
      </tr>
      <tr>
        <td id="L92" class="blob-num js-line-number" data-line-number="92"></td>
        <td id="LC92" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">updateMap</span>(<span class="pl-k">int</span> *cellImage, <span class="pl-k">const</span> <span class="pl-k">int</span> nGlobal, CELLPROP *cellProp, <span class="pl-k">const</span> <span class="pl-k">int</span> nCells);</td>
      </tr>
      <tr>
        <td id="L93" class="blob-num js-line-number" data-line-number="93"></td>
        <td id="LC93" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L94" class="blob-num js-line-number" data-line-number="94"></td>
        <td id="LC94" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L95" class="blob-num js-line-number" data-line-number="95"></td>
        <td id="LC95" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L96" class="blob-num js-line-number" data-line-number="96"></td>
        <td id="LC96" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L97" class="blob-num js-line-number" data-line-number="97"></td>
        <td id="LC97" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L98" class="blob-num js-line-number" data-line-number="98"></td>
        <td id="LC98" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L99" class="blob-num js-line-number" data-line-number="99"></td>
        <td id="LC99" class="blob-code js-file-line"><span class="pl-c">//              vol2bird options from options.conf               //</span></td>
      </tr>
      <tr>
        <td id="L100" class="blob-num js-line-number" data-line-number="100"></td>
        <td id="LC100" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L101" class="blob-num js-line-number" data-line-number="101"></td>
        <td id="LC101" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L102" class="blob-num js-line-number" data-line-number="102"></td>
        <td id="LC102" class="blob-code js-file-line"><span class="pl-c">// the configuration options </span></td>
      </tr>
      <tr>
        <td id="L103" class="blob-num js-line-number" data-line-number="103"></td>
        <td id="LC103" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-c1">cfg_t</span>* cfg;</td>
      </tr>
      <tr>
        <td id="L104" class="blob-num js-line-number" data-line-number="104"></td>
        <td id="LC104" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L105" class="blob-num js-line-number" data-line-number="105"></td>
        <td id="LC105" class="blob-code js-file-line"><span class="pl-c">// the user can specify to exclude gates based on their azimuth;</span></td>
      </tr>
      <tr>
        <td id="L106" class="blob-num js-line-number" data-line-number="106"></td>
        <td id="LC106" class="blob-code js-file-line"><span class="pl-c">// the maximum is set by azimMax</span></td>
      </tr>
      <tr>
        <td id="L107" class="blob-num js-line-number" data-line-number="107"></td>
        <td id="LC107" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> azimMax;</td>
      </tr>
      <tr>
        <td id="L108" class="blob-num js-line-number" data-line-number="108"></td>
        <td id="LC108" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L109" class="blob-num js-line-number" data-line-number="109"></td>
        <td id="LC109" class="blob-code js-file-line"><span class="pl-c">// the user can specify to exclude gates based on their azimuth;</span></td>
      </tr>
      <tr>
        <td id="L110" class="blob-num js-line-number" data-line-number="110"></td>
        <td id="LC110" class="blob-code js-file-line"><span class="pl-c">// the minimum is set by azimMin</span></td>
      </tr>
      <tr>
        <td id="L111" class="blob-num js-line-number" data-line-number="111"></td>
        <td id="LC111" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> azimMin;</td>
      </tr>
      <tr>
        <td id="L112" class="blob-num js-line-number" data-line-number="112"></td>
        <td id="LC112" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L113" class="blob-num js-line-number" data-line-number="113"></td>
        <td id="LC113" class="blob-code js-file-line"><span class="pl-c">// the thickness in meters of a layer in the altitude profile</span></td>
      </tr>
      <tr>
        <td id="L114" class="blob-num js-line-number" data-line-number="114"></td>
        <td id="LC114" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> layerThickness;</td>
      </tr>
      <tr>
        <td id="L115" class="blob-num js-line-number" data-line-number="115"></td>
        <td id="LC115" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L116" class="blob-num js-line-number" data-line-number="116"></td>
        <td id="LC116" class="blob-code js-file-line"><span class="pl-c">// the number of layers in an altitude profile</span></td>
      </tr>
      <tr>
        <td id="L117" class="blob-num js-line-number" data-line-number="117"></td>
        <td id="LC117" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nLayers;</td>
      </tr>
      <tr>
        <td id="L118" class="blob-num js-line-number" data-line-number="118"></td>
        <td id="LC118" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L119" class="blob-num js-line-number" data-line-number="119"></td>
        <td id="LC119" class="blob-code js-file-line"><span class="pl-c">// whether or not to print each scan&#39;s raw dbz values to stderr</span></td>
      </tr>
      <tr>
        <td id="L120" class="blob-num js-line-number" data-line-number="120"></td>
        <td id="LC120" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printDbz;</td>
      </tr>
      <tr>
        <td id="L121" class="blob-num js-line-number" data-line-number="121"></td>
        <td id="LC121" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L122" class="blob-num js-line-number" data-line-number="122"></td>
        <td id="LC122" class="blob-code js-file-line"><span class="pl-c">// whether or not to print each scan&#39;s raw vrad values to stderr</span></td>
      </tr>
      <tr>
        <td id="L123" class="blob-num js-line-number" data-line-number="123"></td>
        <td id="LC123" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printVrad;</td>
      </tr>
      <tr>
        <td id="L124" class="blob-num js-line-number" data-line-number="124"></td>
        <td id="LC124" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L125" class="blob-num js-line-number" data-line-number="125"></td>
        <td id="LC125" class="blob-code js-file-line"><span class="pl-c">// whether or not to print each scan&#39;s raw tex values to stderr</span></td>
      </tr>
      <tr>
        <td id="L126" class="blob-num js-line-number" data-line-number="126"></td>
        <td id="LC126" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printTex;</td>
      </tr>
      <tr>
        <td id="L127" class="blob-num js-line-number" data-line-number="127"></td>
        <td id="LC127" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L128" class="blob-num js-line-number" data-line-number="128"></td>
        <td id="LC128" class="blob-code js-file-line"><span class="pl-c">// whether or not to print each scan&#39;s cell values to stderr</span></td>
      </tr>
      <tr>
        <td id="L129" class="blob-num js-line-number" data-line-number="129"></td>
        <td id="LC129" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printCell;</td>
      </tr>
      <tr>
        <td id="L130" class="blob-num js-line-number" data-line-number="130"></td>
        <td id="LC130" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L131" class="blob-num js-line-number" data-line-number="131"></td>
        <td id="LC131" class="blob-code js-file-line"><span class="pl-c">// whether or not to print each scan&#39;s cell properties to stderr</span></td>
      </tr>
      <tr>
        <td id="L132" class="blob-num js-line-number" data-line-number="132"></td>
        <td id="LC132" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printCellProp;</td>
      </tr>
      <tr>
        <td id="L133" class="blob-num js-line-number" data-line-number="133"></td>
        <td id="LC133" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L134" class="blob-num js-line-number" data-line-number="134"></td>
        <td id="LC134" class="blob-code js-file-line"><span class="pl-c">// whether or not to print each scan&#39;s clutter values to stderr</span></td>
      </tr>
      <tr>
        <td id="L135" class="blob-num js-line-number" data-line-number="135"></td>
        <td id="LC135" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printClut;</td>
      </tr>
      <tr>
        <td id="L136" class="blob-num js-line-number" data-line-number="136"></td>
        <td id="LC136" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L137" class="blob-num js-line-number" data-line-number="137"></td>
        <td id="LC137" class="blob-code js-file-line"><span class="pl-c">// whether or not to print vol2bird&#39;s configuration options to stderr</span></td>
      </tr>
      <tr>
        <td id="L138" class="blob-num js-line-number" data-line-number="138"></td>
        <td id="LC138" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printOptions;</td>
      </tr>
      <tr>
        <td id="L139" class="blob-num js-line-number" data-line-number="139"></td>
        <td id="LC139" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L140" class="blob-num js-line-number" data-line-number="140"></td>
        <td id="LC140" class="blob-code js-file-line"><span class="pl-c">// whether or not to print vol2bird&#39;s profiled data to stderr</span></td>
      </tr>
      <tr>
        <td id="L141" class="blob-num js-line-number" data-line-number="141"></td>
        <td id="LC141" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printProfileVar;</td>
      </tr>
      <tr>
        <td id="L142" class="blob-num js-line-number" data-line-number="142"></td>
        <td id="LC142" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L143" class="blob-num js-line-number" data-line-number="143"></td>
        <td id="LC143" class="blob-code js-file-line"><span class="pl-c">// whether or not to print the &#39;points&#39; array</span></td>
      </tr>
      <tr>
        <td id="L144" class="blob-num js-line-number" data-line-number="144"></td>
        <td id="LC144" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> printPointsArray;</td>
      </tr>
      <tr>
        <td id="L145" class="blob-num js-line-number" data-line-number="145"></td>
        <td id="LC145" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L146" class="blob-num js-line-number" data-line-number="146"></td>
        <td id="LC146" class="blob-code js-file-line"><span class="pl-c">// whether or not a model should be fitted to the observed vrad values</span></td>
      </tr>
      <tr>
        <td id="L147" class="blob-num js-line-number" data-line-number="147"></td>
        <td id="LC147" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> fitVrad;</td>
      </tr>
      <tr>
        <td id="L148" class="blob-num js-line-number" data-line-number="148"></td>
        <td id="LC148" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L149" class="blob-num js-line-number" data-line-number="149"></td>
        <td id="LC149" class="blob-code js-file-line"><span class="pl-c">// the range beyond which observations are excluded when constructing </span></td>
      </tr>
      <tr>
        <td id="L150" class="blob-num js-line-number" data-line-number="150"></td>
        <td id="LC150" class="blob-code js-file-line"><span class="pl-c">// the altitude profile</span></td>
      </tr>
      <tr>
        <td id="L151" class="blob-num js-line-number" data-line-number="151"></td>
        <td id="LC151" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> rangeMax;</td>
      </tr>
      <tr>
        <td id="L152" class="blob-num js-line-number" data-line-number="152"></td>
        <td id="LC152" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L153" class="blob-num js-line-number" data-line-number="153"></td>
        <td id="LC153" class="blob-code js-file-line"><span class="pl-c">// the range below which observations are excluded when constructing </span></td>
      </tr>
      <tr>
        <td id="L154" class="blob-num js-line-number" data-line-number="154"></td>
        <td id="LC154" class="blob-code js-file-line"><span class="pl-c">// the altitude profile</span></td>
      </tr>
      <tr>
        <td id="L155" class="blob-num js-line-number" data-line-number="155"></td>
        <td id="LC155" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> rangeMin;</td>
      </tr>
      <tr>
        <td id="L156" class="blob-num js-line-number" data-line-number="156"></td>
        <td id="LC156" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L157" class="blob-num js-line-number" data-line-number="157"></td>
        <td id="LC157" class="blob-code js-file-line"><span class="pl-c">// the wavelength of the radar in units of centimeter</span></td>
      </tr>
      <tr>
        <td id="L158" class="blob-num js-line-number" data-line-number="158"></td>
        <td id="LC158" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> radarWavelength;</td>
      </tr>
      <tr>
        <td id="L159" class="blob-num js-line-number" data-line-number="159"></td>
        <td id="LC159" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L160" class="blob-num js-line-number" data-line-number="160"></td>
        <td id="LC160" class="blob-code js-file-line"><span class="pl-c">// whether clutter data is used</span></td>
      </tr>
      <tr>
        <td id="L161" class="blob-num js-line-number" data-line-number="161"></td>
        <td id="LC161" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> useStaticClutterData;</td>
      </tr>
      <tr>
        <td id="L162" class="blob-num js-line-number" data-line-number="162"></td>
        <td id="LC162" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L163" class="blob-num js-line-number" data-line-number="163"></td>
        <td id="LC163" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L164" class="blob-num js-line-number" data-line-number="164"></td>
        <td id="LC164" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L165" class="blob-num js-line-number" data-line-number="165"></td>
        <td id="LC165" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L166" class="blob-num js-line-number" data-line-number="166"></td>
        <td id="LC166" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L167" class="blob-num js-line-number" data-line-number="167"></td>
        <td id="LC167" class="blob-code js-file-line"><span class="pl-c">//              vol2bird options from constants.h                //</span></td>
      </tr>
      <tr>
        <td id="L168" class="blob-num js-line-number" data-line-number="168"></td>
        <td id="LC168" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L169" class="blob-num js-line-number" data-line-number="169"></td>
        <td id="LC169" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L170" class="blob-num js-line-number" data-line-number="170"></td>
        <td id="LC170" class="blob-code js-file-line"><span class="pl-c">// when analyzing cells, AREAMIN determines the minimum size of a </span></td>
      </tr>
      <tr>
        <td id="L171" class="blob-num js-line-number" data-line-number="171"></td>
        <td id="LC171" class="blob-code js-file-line"><span class="pl-c">// cell to be considered in the rest of the analysis</span></td>
      </tr>
      <tr>
        <td id="L172" class="blob-num js-line-number" data-line-number="172"></td>
        <td id="LC172" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nGatesCellMin;</td>
      </tr>
      <tr>
        <td id="L173" class="blob-num js-line-number" data-line-number="173"></td>
        <td id="LC173" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L174" class="blob-num js-line-number" data-line-number="174"></td>
        <td id="LC174" class="blob-code js-file-line"><span class="pl-c">// Maximum percentage of clutter in a cell </span></td>
      </tr>
      <tr>
        <td id="L175" class="blob-num js-line-number" data-line-number="175"></td>
        <td id="LC175" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> cellClutterFractionMax;</td>
      </tr>
      <tr>
        <td id="L176" class="blob-num js-line-number" data-line-number="176"></td>
        <td id="LC176" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L177" class="blob-num js-line-number" data-line-number="177"></td>
        <td id="LC177" class="blob-code js-file-line"><span class="pl-c">// when analyzing cells, only cells for which the average dbz is </span></td>
      </tr>
      <tr>
        <td id="L178" class="blob-num js-line-number" data-line-number="178"></td>
        <td id="LC178" class="blob-code js-file-line"><span class="pl-c">// more than DBZCELL are considered in the rest of the analysis</span></td>
      </tr>
      <tr>
        <td id="L179" class="blob-num js-line-number" data-line-number="179"></td>
        <td id="LC179" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> cellDbzMin;</td>
      </tr>
      <tr>
        <td id="L180" class="blob-num js-line-number" data-line-number="180"></td>
        <td id="LC180" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L181" class="blob-num js-line-number" data-line-number="181"></td>
        <td id="LC181" class="blob-code js-file-line"><span class="pl-c">// minimum quality of the fit</span></td>
      </tr>
      <tr>
        <td id="L182" class="blob-num js-line-number" data-line-number="182"></td>
        <td id="LC182" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> chisqMin;</td>
      </tr>
      <tr>
        <td id="L183" class="blob-num js-line-number" data-line-number="183"></td>
        <td id="LC183" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L184" class="blob-num js-line-number" data-line-number="184"></td>
        <td id="LC184" class="blob-code js-file-line"><span class="pl-c">// threshold dbz value for excluding gates as clutter (static clutter </span></td>
      </tr>
      <tr>
        <td id="L185" class="blob-num js-line-number" data-line-number="185"></td>
        <td id="LC185" class="blob-code js-file-line"><span class="pl-c">// only)</span></td>
      </tr>
      <tr>
        <td id="L186" class="blob-num js-line-number" data-line-number="186"></td>
        <td id="LC186" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> clutterValueMin;</td>
      </tr>
      <tr>
        <td id="L187" class="blob-num js-line-number" data-line-number="187"></td>
        <td id="LC187" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L188" class="blob-num js-line-number" data-line-number="188"></td>
        <td id="LC188" class="blob-code js-file-line"><span class="pl-c">// maximum dbz used in calculation of profile dbzAvg  </span></td>
      </tr>
      <tr>
        <td id="L189" class="blob-num js-line-number" data-line-number="189"></td>
        <td id="LC189" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> dbzMax;</td>
      </tr>
      <tr>
        <td id="L190" class="blob-num js-line-number" data-line-number="190"></td>
        <td id="LC190" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L191" class="blob-num js-line-number" data-line-number="191"></td>
        <td id="LC191" class="blob-code js-file-line"><span class="pl-c">// minimum dbz used in calculation of cell dbzAvg</span></td>
      </tr>
      <tr>
        <td id="L192" class="blob-num js-line-number" data-line-number="192"></td>
        <td id="LC192" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> dbzThresMin;</td>
      </tr>
      <tr>
        <td id="L193" class="blob-num js-line-number" data-line-number="193"></td>
        <td id="LC193" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L194" class="blob-num js-line-number" data-line-number="194"></td>
        <td id="LC194" class="blob-code js-file-line"><span class="pl-c">// each weather cell identified by findWeatherCells() is grown by a distance </span></td>
      </tr>
      <tr>
        <td id="L195" class="blob-num js-line-number" data-line-number="195"></td>
        <td id="LC195" class="blob-code js-file-line"><span class="pl-c">// equal to &#39;fringeDist&#39; using a region-growing approach</span></td>
      </tr>
      <tr>
        <td id="L196" class="blob-num js-line-number" data-line-number="196"></td>
        <td id="LC196" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> fringeDist;</td>
      </tr>
      <tr>
        <td id="L197" class="blob-num js-line-number" data-line-number="197"></td>
        <td id="LC197" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L198" class="blob-num js-line-number" data-line-number="198"></td>
        <td id="LC198" class="blob-code js-file-line"><span class="pl-c">// when determining whether there are enough vrad observations in </span></td>
      </tr>
      <tr>
        <td id="L199" class="blob-num js-line-number" data-line-number="199"></td>
        <td id="LC199" class="blob-code js-file-line"><span class="pl-c">// each direction, use NBINSGAP sectors</span></td>
      </tr>
      <tr>
        <td id="L200" class="blob-num js-line-number" data-line-number="200"></td>
        <td id="LC200" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nBinsGap;</td>
      </tr>
      <tr>
        <td id="L201" class="blob-num js-line-number" data-line-number="201"></td>
        <td id="LC201" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L202" class="blob-num js-line-number" data-line-number="202"></td>
        <td id="LC202" class="blob-code js-file-line"><span class="pl-c">// when calculating the altitude-layer averaged dbz, there should </span></td>
      </tr>
      <tr>
        <td id="L203" class="blob-num js-line-number" data-line-number="203"></td>
        <td id="LC203" class="blob-code js-file-line"><span class="pl-c">// be at least NDBZMIN valid data points</span></td>
      </tr>
      <tr>
        <td id="L204" class="blob-num js-line-number" data-line-number="204"></td>
        <td id="LC204" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nPointsIncludedMin;</td>
      </tr>
      <tr>
        <td id="L205" class="blob-num js-line-number" data-line-number="205"></td>
        <td id="LC205" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L206" class="blob-num js-line-number" data-line-number="206"></td>
        <td id="LC206" class="blob-code js-file-line"><span class="pl-c">// the minimum number of direct neighbors with dbz value above </span></td>
      </tr>
      <tr>
        <td id="L207" class="blob-num js-line-number" data-line-number="207"></td>
        <td id="LC207" class="blob-code js-file-line"><span class="pl-c">// dbzThresMin as used in findWeatherCells()</span></td>
      </tr>
      <tr>
        <td id="L208" class="blob-num js-line-number" data-line-number="208"></td>
        <td id="LC208" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nNeighborsMin;</td>
      </tr>
      <tr>
        <td id="L209" class="blob-num js-line-number" data-line-number="209"></td>
        <td id="LC209" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L210" class="blob-num js-line-number" data-line-number="210"></td>
        <td id="LC210" class="blob-code js-file-line"><span class="pl-c">// there should be at least NOBSGAPMIN vrad observations in each </span></td>
      </tr>
      <tr>
        <td id="L211" class="blob-num js-line-number" data-line-number="211"></td>
        <td id="LC211" class="blob-code js-file-line"><span class="pl-c">// sector</span></td>
      </tr>
      <tr>
        <td id="L212" class="blob-num js-line-number" data-line-number="212"></td>
        <td id="LC212" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nObsGapMin;</td>
      </tr>
      <tr>
        <td id="L213" class="blob-num js-line-number" data-line-number="213"></td>
        <td id="LC213" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L214" class="blob-num js-line-number" data-line-number="214"></td>
        <td id="LC214" class="blob-code js-file-line"><span class="pl-c">// vrad&#39;s texture is calculated based on the local neighborhood. The</span></td>
      </tr>
      <tr>
        <td id="L215" class="blob-num js-line-number" data-line-number="215"></td>
        <td id="LC215" class="blob-code js-file-line"><span class="pl-c">// neighborhood size in the azimuth direction is equal to NTEXBINAZIM</span></td>
      </tr>
      <tr>
        <td id="L216" class="blob-num js-line-number" data-line-number="216"></td>
        <td id="LC216" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nAzimNeighborhood;</td>
      </tr>
      <tr>
        <td id="L217" class="blob-num js-line-number" data-line-number="217"></td>
        <td id="LC217" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L218" class="blob-num js-line-number" data-line-number="218"></td>
        <td id="LC218" class="blob-code js-file-line"><span class="pl-c">// vrad&#39;s texture is calculated based on the local neighborhood. The</span></td>
      </tr>
      <tr>
        <td id="L219" class="blob-num js-line-number" data-line-number="219"></td>
        <td id="LC219" class="blob-code js-file-line"><span class="pl-c">// neighborhood size in the range direction is equal to NTEXBINRANG</span></td>
      </tr>
      <tr>
        <td id="L220" class="blob-num js-line-number" data-line-number="220"></td>
        <td id="LC220" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nRangNeighborhood;</td>
      </tr>
      <tr>
        <td id="L221" class="blob-num js-line-number" data-line-number="221"></td>
        <td id="LC221" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L222" class="blob-num js-line-number" data-line-number="222"></td>
        <td id="LC222" class="blob-code js-file-line"><span class="pl-c">// the minimum number of neighbors for the texture value to be </span></td>
      </tr>
      <tr>
        <td id="L223" class="blob-num js-line-number" data-line-number="223"></td>
        <td id="LC223" class="blob-code js-file-line"><span class="pl-c">// considered valid, as used in calcTexture()</span></td>
      </tr>
      <tr>
        <td id="L224" class="blob-num js-line-number" data-line-number="224"></td>
        <td id="LC224" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nCountMin; </td>
      </tr>
      <tr>
        <td id="L225" class="blob-num js-line-number" data-line-number="225"></td>
        <td id="LC225" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L226" class="blob-num js-line-number" data-line-number="226"></td>
        <td id="LC226" class="blob-code js-file-line"><span class="pl-c">// the refractive index of water</span></td>
      </tr>
      <tr>
        <td id="L227" class="blob-num js-line-number" data-line-number="227"></td>
        <td id="LC227" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> refracIndex;</td>
      </tr>
      <tr>
        <td id="L228" class="blob-num js-line-number" data-line-number="228"></td>
        <td id="LC228" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L229" class="blob-num js-line-number" data-line-number="229"></td>
        <td id="LC229" class="blob-code js-file-line"><span class="pl-c">// the bird radar cross section</span></td>
      </tr>
      <tr>
        <td id="L230" class="blob-num js-line-number" data-line-number="230"></td>
        <td id="LC230" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> birdRadarCrossSection;</td>
      </tr>
      <tr>
        <td id="L231" class="blob-num js-line-number" data-line-number="231"></td>
        <td id="LC231" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L232" class="blob-num js-line-number" data-line-number="232"></td>
        <td id="LC232" class="blob-code js-file-line"><span class="pl-c">// when analyzing cells, only cells for which the stddev of vrad</span></td>
      </tr>
      <tr>
        <td id="L233" class="blob-num js-line-number" data-line-number="233"></td>
        <td id="LC233" class="blob-code js-file-line"><span class="pl-c">// (aka the texture) is less than STDEVCELL are considered in the</span></td>
      </tr>
      <tr>
        <td id="L234" class="blob-num js-line-number" data-line-number="234"></td>
        <td id="LC234" class="blob-code js-file-line"><span class="pl-c">// rest of the analysis</span></td>
      </tr>
      <tr>
        <td id="L235" class="blob-num js-line-number" data-line-number="235"></td>
        <td id="LC235" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> cellStdDevMax;</td>
      </tr>
      <tr>
        <td id="L236" class="blob-num js-line-number" data-line-number="236"></td>
        <td id="LC236" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L237" class="blob-num js-line-number" data-line-number="237"></td>
        <td id="LC237" class="blob-code js-file-line"><span class="pl-c">// Minimum radial velocity standard deviation to be labeled as birds</span></td>
      </tr>
      <tr>
        <td id="L238" class="blob-num js-line-number" data-line-number="238"></td>
        <td id="LC238" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> stdDevMinBird;</td>
      </tr>
      <tr>
        <td id="L239" class="blob-num js-line-number" data-line-number="239"></td>
        <td id="LC239" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L240" class="blob-num js-line-number" data-line-number="240"></td>
        <td id="LC240" class="blob-code js-file-line"><span class="pl-c">// after fitting the vrad data, throw out any vrad observations that </span></td>
      </tr>
      <tr>
        <td id="L241" class="blob-num js-line-number" data-line-number="241"></td>
        <td id="LC241" class="blob-code js-file-line"><span class="pl-c">// are more that VDIFMAX away from the fitted value, since these are</span></td>
      </tr>
      <tr>
        <td id="L242" class="blob-num js-line-number" data-line-number="242"></td>
        <td id="LC242" class="blob-code js-file-line"><span class="pl-c">// likely outliers </span></td>
      </tr>
      <tr>
        <td id="L243" class="blob-num js-line-number" data-line-number="243"></td>
        <td id="LC243" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> absVDifMax;</td>
      </tr>
      <tr>
        <td id="L244" class="blob-num js-line-number" data-line-number="244"></td>
        <td id="LC244" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L245" class="blob-num js-line-number" data-line-number="245"></td>
        <td id="LC245" class="blob-code js-file-line"><span class="pl-c">// When analyzing cells, radial velocities lower than VRADMIN </span></td>
      </tr>
      <tr>
        <td id="L246" class="blob-num js-line-number" data-line-number="246"></td>
        <td id="LC246" class="blob-code js-file-line"><span class="pl-c">// are treated as clutter</span></td>
      </tr>
      <tr>
        <td id="L247" class="blob-num js-line-number" data-line-number="247"></td>
        <td id="LC247" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> vradMin;</td>
      </tr>
      <tr>
        <td id="L248" class="blob-num js-line-number" data-line-number="248"></td>
        <td id="LC248" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L249" class="blob-num js-line-number" data-line-number="249"></td>
        <td id="LC249" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L250" class="blob-num js-line-number" data-line-number="250"></td>
        <td id="LC250" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L251" class="blob-num js-line-number" data-line-number="251"></td>
        <td id="LC251" class="blob-code js-file-line"><span class="pl-c">//               information about the &#39;points&#39; array            //</span></td>
      </tr>
      <tr>
        <td id="L252" class="blob-num js-line-number" data-line-number="252"></td>
        <td id="LC252" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L253" class="blob-num js-line-number" data-line-number="253"></td>
        <td id="LC253" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L254" class="blob-num js-line-number" data-line-number="254"></td>
        <td id="LC254" class="blob-code js-file-line"><span class="pl-c">// The data needed for calculating bird densities are collected</span></td>
      </tr>
      <tr>
        <td id="L255" class="blob-num js-line-number" data-line-number="255"></td>
        <td id="LC255" class="blob-code js-file-line"><span class="pl-c">// in one big array, &#39;points&#39;. Although this array is one </span></td>
      </tr>
      <tr>
        <td id="L256" class="blob-num js-line-number" data-line-number="256"></td>
        <td id="LC256" class="blob-code js-file-line"><span class="pl-c">// variable, it is partitioned into &#39;nLayers&#39; parts. The parts </span></td>
      </tr>
      <tr>
        <td id="L257" class="blob-num js-line-number" data-line-number="257"></td>
        <td id="LC257" class="blob-code js-file-line"><span class="pl-c">// are not equal in size, therefore we need to keep track of  </span></td>
      </tr>
      <tr>
        <td id="L258" class="blob-num js-line-number" data-line-number="258"></td>
        <td id="LC258" class="blob-code js-file-line"><span class="pl-c">// where the data pertaining to a certain altitude bin can be</span></td>
      </tr>
      <tr>
        <td id="L259" class="blob-num js-line-number" data-line-number="259"></td>
        <td id="LC259" class="blob-code js-file-line"><span class="pl-c">// written. The valid range of indexes into &#39;points&#39; are stored </span></td>
      </tr>
      <tr>
        <td id="L260" class="blob-num js-line-number" data-line-number="260"></td>
        <td id="LC260" class="blob-code js-file-line"><span class="pl-c">// in arrays &#39;indexFrom&#39; and &#39;indexTo&#39;.</span></td>
      </tr>
      <tr>
        <td id="L261" class="blob-num js-line-number" data-line-number="261"></td>
        <td id="LC261" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L262" class="blob-num js-line-number" data-line-number="262"></td>
        <td id="LC262" class="blob-code js-file-line"><span class="pl-c">// the &#39;points&#39; array has this many pseudo-columns</span></td>
      </tr>
      <tr>
        <td id="L263" class="blob-num js-line-number" data-line-number="263"></td>
        <td id="LC263" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nColsPoints;</td>
      </tr>
      <tr>
        <td id="L264" class="blob-num js-line-number" data-line-number="264"></td>
        <td id="LC264" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L265" class="blob-num js-line-number" data-line-number="265"></td>
        <td id="LC265" class="blob-code js-file-line"><span class="pl-c">// the &#39;points&#39; array has this many rows</span></td>
      </tr>
      <tr>
        <td id="L266" class="blob-num js-line-number" data-line-number="266"></td>
        <td id="LC266" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nRowsPoints;</td>
      </tr>
      <tr>
        <td id="L267" class="blob-num js-line-number" data-line-number="267"></td>
        <td id="LC267" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L268" class="blob-num js-line-number" data-line-number="268"></td>
        <td id="LC268" class="blob-code js-file-line"><span class="pl-c">// the psuedo-column in &#39;points&#39; that holds the azimuth angle</span></td>
      </tr>
      <tr>
        <td id="L269" class="blob-num js-line-number" data-line-number="269"></td>
        <td id="LC269" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> azimAngleCol;</td>
      </tr>
      <tr>
        <td id="L270" class="blob-num js-line-number" data-line-number="270"></td>
        <td id="LC270" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L271" class="blob-num js-line-number" data-line-number="271"></td>
        <td id="LC271" class="blob-code js-file-line"><span class="pl-c">// the psuedo-column in &#39;points&#39; that holds the elevation angle</span></td>
      </tr>
      <tr>
        <td id="L272" class="blob-num js-line-number" data-line-number="272"></td>
        <td id="LC272" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> elevAngleCol;</td>
      </tr>
      <tr>
        <td id="L273" class="blob-num js-line-number" data-line-number="273"></td>
        <td id="LC273" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L274" class="blob-num js-line-number" data-line-number="274"></td>
        <td id="LC274" class="blob-code js-file-line"><span class="pl-c">// the psuedo-column in &#39;points&#39; that holds the dbz value</span></td>
      </tr>
      <tr>
        <td id="L275" class="blob-num js-line-number" data-line-number="275"></td>
        <td id="LC275" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> dbzValueCol;</td>
      </tr>
      <tr>
        <td id="L276" class="blob-num js-line-number" data-line-number="276"></td>
        <td id="LC276" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L277" class="blob-num js-line-number" data-line-number="277"></td>
        <td id="LC277" class="blob-code js-file-line"><span class="pl-c">// the psuedo-column in &#39;points&#39; that holds the vrad value</span></td>
      </tr>
      <tr>
        <td id="L278" class="blob-num js-line-number" data-line-number="278"></td>
        <td id="LC278" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> vradValueCol;</td>
      </tr>
      <tr>
        <td id="L279" class="blob-num js-line-number" data-line-number="279"></td>
        <td id="LC279" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L280" class="blob-num js-line-number" data-line-number="280"></td>
        <td id="LC280" class="blob-code js-file-line"><span class="pl-c">// the psuedo-column in &#39;points&#39; that holds the cell value</span></td>
      </tr>
      <tr>
        <td id="L281" class="blob-num js-line-number" data-line-number="281"></td>
        <td id="LC281" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> cellValueCol;</td>
      </tr>
      <tr>
        <td id="L282" class="blob-num js-line-number" data-line-number="282"></td>
        <td id="LC282" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L283" class="blob-num js-line-number" data-line-number="283"></td>
        <td id="LC283" class="blob-code js-file-line"><span class="pl-c">// the psuedo-column in &#39;points&#39; that holds the gate classification code</span></td>
      </tr>
      <tr>
        <td id="L284" class="blob-num js-line-number" data-line-number="284"></td>
        <td id="LC284" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> gateCodeCol;</td>
      </tr>
      <tr>
        <td id="L285" class="blob-num js-line-number" data-line-number="285"></td>
        <td id="LC285" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L286" class="blob-num js-line-number" data-line-number="286"></td>
        <td id="LC286" class="blob-code js-file-line"><span class="pl-c">// the &#39;points&#39; array itself</span></td>
      </tr>
      <tr>
        <td id="L287" class="blob-num js-line-number" data-line-number="287"></td>
        <td id="LC287" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span>* points;</td>
      </tr>
      <tr>
        <td id="L288" class="blob-num js-line-number" data-line-number="288"></td>
        <td id="LC288" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L289" class="blob-num js-line-number" data-line-number="289"></td>
        <td id="LC289" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L290" class="blob-num js-line-number" data-line-number="290"></td>
        <td id="LC290" class="blob-code js-file-line"><span class="pl-c">//             lists of indices into the &#39;points&#39; array:         //</span></td>
      </tr>
      <tr>
        <td id="L291" class="blob-num js-line-number" data-line-number="291"></td>
        <td id="LC291" class="blob-code js-file-line"><span class="pl-c">//          where each altitude layer&#39;s data starts and ends     //</span></td>
      </tr>
      <tr>
        <td id="L292" class="blob-num js-line-number" data-line-number="292"></td>
        <td id="LC292" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L293" class="blob-num js-line-number" data-line-number="293"></td>
        <td id="LC293" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L294" class="blob-num js-line-number" data-line-number="294"></td>
        <td id="LC294" class="blob-code js-file-line"><span class="pl-c">// for a given altitude layer in the profile, only part of the &#39;points&#39;</span></td>
      </tr>
      <tr>
        <td id="L295" class="blob-num js-line-number" data-line-number="295"></td>
        <td id="LC295" class="blob-code js-file-line"><span class="pl-c">// array is relevant. The &#39;indexFrom&#39; and &#39;indexTo&#39; arrays keep track </span></td>
      </tr>
      <tr>
        <td id="L296" class="blob-num js-line-number" data-line-number="296"></td>
        <td id="LC296" class="blob-code js-file-line"><span class="pl-c">// which rows in &#39;points&#39; pertains to a given layer</span></td>
      </tr>
      <tr>
        <td id="L297" class="blob-num js-line-number" data-line-number="297"></td>
        <td id="LC297" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span>* indexFrom;</td>
      </tr>
      <tr>
        <td id="L298" class="blob-num js-line-number" data-line-number="298"></td>
        <td id="LC298" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span>* indexTo;</td>
      </tr>
      <tr>
        <td id="L299" class="blob-num js-line-number" data-line-number="299"></td>
        <td id="LC299" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L300" class="blob-num js-line-number" data-line-number="300"></td>
        <td id="LC300" class="blob-code js-file-line"><span class="pl-c">// nPointsWritten stores the number of points that was copied from one </span></td>
      </tr>
      <tr>
        <td id="L301" class="blob-num js-line-number" data-line-number="301"></td>
        <td id="LC301" class="blob-code js-file-line"><span class="pl-c">// of the scan elevations to the &#39;points&#39; array; it should therefore</span></td>
      </tr>
      <tr>
        <td id="L302" class="blob-num js-line-number" data-line-number="302"></td>
        <td id="LC302" class="blob-code js-file-line"><span class="pl-c">// never exceed indexTo[i]-indexFrom[i]</span></td>
      </tr>
      <tr>
        <td id="L303" class="blob-num js-line-number" data-line-number="303"></td>
        <td id="LC303" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span>* nPointsWritten;</td>
      </tr>
      <tr>
        <td id="L304" class="blob-num js-line-number" data-line-number="304"></td>
        <td id="LC304" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L305" class="blob-num js-line-number" data-line-number="305"></td>
        <td id="LC305" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L306" class="blob-num js-line-number" data-line-number="306"></td>
        <td id="LC306" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L307" class="blob-num js-line-number" data-line-number="307"></td>
        <td id="LC307" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L308" class="blob-num js-line-number" data-line-number="308"></td>
        <td id="LC308" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L309" class="blob-num js-line-number" data-line-number="309"></td>
        <td id="LC309" class="blob-code js-file-line"><span class="pl-c">//          information about the flagfields of &#39;gateCode&#39;       //</span></td>
      </tr>
      <tr>
        <td id="L310" class="blob-num js-line-number" data-line-number="310"></td>
        <td id="LC310" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L311" class="blob-num js-line-number" data-line-number="311"></td>
        <td id="LC311" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L312" class="blob-num js-line-number" data-line-number="312"></td>
        <td id="LC312" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; that says whether this gate is true in the static</span></td>
      </tr>
      <tr>
        <td id="L313" class="blob-num js-line-number" data-line-number="313"></td>
        <td id="LC313" class="blob-code js-file-line"><span class="pl-c">// clutter map (which we don&#39;t have yet TODO)</span></td>
      </tr>
      <tr>
        <td id="L314" class="blob-num js-line-number" data-line-number="314"></td>
        <td id="LC314" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionStaticClutter;</td>
      </tr>
      <tr>
        <td id="L315" class="blob-num js-line-number" data-line-number="315"></td>
        <td id="LC315" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L316" class="blob-num js-line-number" data-line-number="316"></td>
        <td id="LC316" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; that says whether this gate is part of the </span></td>
      </tr>
      <tr>
        <td id="L317" class="blob-num js-line-number" data-line-number="317"></td>
        <td id="LC317" class="blob-code js-file-line"><span class="pl-c">// calculated cluttermap (without fringe)</span></td>
      </tr>
      <tr>
        <td id="L318" class="blob-num js-line-number" data-line-number="318"></td>
        <td id="LC318" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionDynamicClutter;</td>
      </tr>
      <tr>
        <td id="L319" class="blob-num js-line-number" data-line-number="319"></td>
        <td id="LC319" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L320" class="blob-num js-line-number" data-line-number="320"></td>
        <td id="LC320" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; that says whether this gate is part of the </span></td>
      </tr>
      <tr>
        <td id="L321" class="blob-num js-line-number" data-line-number="321"></td>
        <td id="LC321" class="blob-code js-file-line"><span class="pl-c">// fringe of the calculated cluttermap</span></td>
      </tr>
      <tr>
        <td id="L322" class="blob-num js-line-number" data-line-number="322"></td>
        <td id="LC322" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionDynamicClutterFringe;</td>
      </tr>
      <tr>
        <td id="L323" class="blob-num js-line-number" data-line-number="323"></td>
        <td id="LC323" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L324" class="blob-num js-line-number" data-line-number="324"></td>
        <td id="LC324" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; that says whether this gate has reflectivity data </span></td>
      </tr>
      <tr>
        <td id="L325" class="blob-num js-line-number" data-line-number="325"></td>
        <td id="LC325" class="blob-code js-file-line"><span class="pl-c">// but no corresponding radial velocity data</span></td>
      </tr>
      <tr>
        <td id="L326" class="blob-num js-line-number" data-line-number="326"></td>
        <td id="LC326" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionVradMissing;</td>
      </tr>
      <tr>
        <td id="L327" class="blob-num js-line-number" data-line-number="327"></td>
        <td id="LC327" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L328" class="blob-num js-line-number" data-line-number="328"></td>
        <td id="LC328" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; the psuedo-columnsays whether this gate&#39;s dbz value is too</span></td>
      </tr>
      <tr>
        <td id="L329" class="blob-num js-line-number" data-line-number="329"></td>
        <td id="LC329" class="blob-code js-file-line"><span class="pl-c">// high to be due to birds, it must be caused by something else</span></td>
      </tr>
      <tr>
        <td id="L330" class="blob-num js-line-number" data-line-number="330"></td>
        <td id="LC330" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionDbzTooHighForBirds;</td>
      </tr>
      <tr>
        <td id="L331" class="blob-num js-line-number" data-line-number="331"></td>
        <td id="LC331" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L332" class="blob-num js-line-number" data-line-number="332"></td>
        <td id="LC332" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; the psuedo-columnsays whether this gate&#39;s radial velocity is</span></td>
      </tr>
      <tr>
        <td id="L333" class="blob-num js-line-number" data-line-number="333"></td>
        <td id="LC333" class="blob-code js-file-line"><span class="pl-c">// close to zero. These gates are all discarded to exclude ground </span></td>
      </tr>
      <tr>
        <td id="L334" class="blob-num js-line-number" data-line-number="334"></td>
        <td id="LC334" class="blob-code js-file-line"><span class="pl-c">// clutter, which often has a radial velocity near zero.</span></td>
      </tr>
      <tr>
        <td id="L335" class="blob-num js-line-number" data-line-number="335"></td>
        <td id="LC335" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionVradTooLow;</td>
      </tr>
      <tr>
        <td id="L336" class="blob-num js-line-number" data-line-number="336"></td>
        <td id="LC336" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L337" class="blob-num js-line-number" data-line-number="337"></td>
        <td id="LC337" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; that says whether this gate passed the VDIFMAX test</span></td>
      </tr>
      <tr>
        <td id="L338" class="blob-num js-line-number" data-line-number="338"></td>
        <td id="LC338" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionVDifMax;</td>
      </tr>
      <tr>
        <td id="L339" class="blob-num js-line-number" data-line-number="339"></td>
        <td id="LC339" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L340" class="blob-num js-line-number" data-line-number="340"></td>
        <td id="LC340" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; that says whether the gate&#39;s azimuth angle was too low</span></td>
      </tr>
      <tr>
        <td id="L341" class="blob-num js-line-number" data-line-number="341"></td>
        <td id="LC341" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionAzimTooLow;</td>
      </tr>
      <tr>
        <td id="L342" class="blob-num js-line-number" data-line-number="342"></td>
        <td id="LC342" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L343" class="blob-num js-line-number" data-line-number="343"></td>
        <td id="LC343" class="blob-code js-file-line"><span class="pl-c">// the bit in &#39;gateCode&#39; that says whether the gate&#39;s azimuth angle was too high</span></td>
      </tr>
      <tr>
        <td id="L344" class="blob-num js-line-number" data-line-number="344"></td>
        <td id="LC344" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> flagPositionAzimTooHigh;</td>
      </tr>
      <tr>
        <td id="L345" class="blob-num js-line-number" data-line-number="345"></td>
        <td id="LC345" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L346" class="blob-num js-line-number" data-line-number="346"></td>
        <td id="LC346" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L347" class="blob-num js-line-number" data-line-number="347"></td>
        <td id="LC347" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L348" class="blob-num js-line-number" data-line-number="348"></td>
        <td id="LC348" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L349" class="blob-num js-line-number" data-line-number="349"></td>
        <td id="LC349" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L350" class="blob-num js-line-number" data-line-number="350"></td>
        <td id="LC350" class="blob-code js-file-line"><span class="pl-c">//              information about the &#39;profile&#39; array            //</span></td>
      </tr>
      <tr>
        <td id="L351" class="blob-num js-line-number" data-line-number="351"></td>
        <td id="LC351" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L352" class="blob-num js-line-number" data-line-number="352"></td>
        <td id="LC352" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L353" class="blob-num js-line-number" data-line-number="353"></td>
        <td id="LC353" class="blob-code js-file-line"><span class="pl-c">// the number of different types of profile we&#39;re making</span></td>
      </tr>
      <tr>
        <td id="L354" class="blob-num js-line-number" data-line-number="354"></td>
        <td id="LC354" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nProfileTypes;</td>
      </tr>
      <tr>
        <td id="L355" class="blob-num js-line-number" data-line-number="355"></td>
        <td id="LC355" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L356" class="blob-num js-line-number" data-line-number="356"></td>
        <td id="LC356" class="blob-code js-file-line"><span class="pl-c">// how many rows there are in a profile</span></td>
      </tr>
      <tr>
        <td id="L357" class="blob-num js-line-number" data-line-number="357"></td>
        <td id="LC357" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nRowsProfile;</td>
      </tr>
      <tr>
        <td id="L358" class="blob-num js-line-number" data-line-number="358"></td>
        <td id="LC358" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L359" class="blob-num js-line-number" data-line-number="359"></td>
        <td id="LC359" class="blob-code js-file-line"><span class="pl-c">// columns in profile contain [altmin,altmax,u,v,w,hSpeed,hDir,chi,hasGap,dbzAvg,nPointsCopied,reflectivity,birdDensity]</span></td>
      </tr>
      <tr>
        <td id="L360" class="blob-num js-line-number" data-line-number="360"></td>
        <td id="LC360" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nColsProfile; </td>
      </tr>
      <tr>
        <td id="L361" class="blob-num js-line-number" data-line-number="361"></td>
        <td id="LC361" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L362" class="blob-num js-line-number" data-line-number="362"></td>
        <td id="LC362" class="blob-code js-file-line"><span class="pl-c">// the profile array itself</span></td>
      </tr>
      <tr>
        <td id="L363" class="blob-num js-line-number" data-line-number="363"></td>
        <td id="LC363" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span>* profile;</td>
      </tr>
      <tr>
        <td id="L364" class="blob-num js-line-number" data-line-number="364"></td>
        <td id="LC364" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L365" class="blob-num js-line-number" data-line-number="365"></td>
        <td id="LC365" class="blob-code js-file-line"><span class="pl-c">// the type of profile that was last calculated</span></td>
      </tr>
      <tr>
        <td id="L366" class="blob-num js-line-number" data-line-number="366"></td>
        <td id="LC366" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> iProfileTypeLast;</td>
      </tr>
      <tr>
        <td id="L367" class="blob-num js-line-number" data-line-number="367"></td>
        <td id="LC367" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L368" class="blob-num js-line-number" data-line-number="368"></td>
        <td id="LC368" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L369" class="blob-num js-line-number" data-line-number="369"></td>
        <td id="LC369" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L370" class="blob-num js-line-number" data-line-number="370"></td>
        <td id="LC370" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L371" class="blob-num js-line-number" data-line-number="371"></td>
        <td id="LC371" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L372" class="blob-num js-line-number" data-line-number="372"></td>
        <td id="LC372" class="blob-code js-file-line"><span class="pl-c">//                       some other variables                    //</span></td>
      </tr>
      <tr>
        <td id="L373" class="blob-num js-line-number" data-line-number="373"></td>
        <td id="LC373" class="blob-code js-file-line"><span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L374" class="blob-num js-line-number" data-line-number="374"></td>
        <td id="LC374" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L375" class="blob-num js-line-number" data-line-number="375"></td>
        <td id="LC375" class="blob-code js-file-line"><span class="pl-c">// rCellMax is defined as rangeMax + 5000.0f  in order to avoid </span></td>
      </tr>
      <tr>
        <td id="L376" class="blob-num js-line-number" data-line-number="376"></td>
        <td id="LC376" class="blob-code js-file-line"><span class="pl-c">// edge effects when calculating the fringe</span></td>
      </tr>
      <tr>
        <td id="L377" class="blob-num js-line-number" data-line-number="377"></td>
        <td id="LC377" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> rCellMax;</td>
      </tr>
      <tr>
        <td id="L378" class="blob-num js-line-number" data-line-number="378"></td>
        <td id="LC378" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L379" class="blob-num js-line-number" data-line-number="379"></td>
        <td id="LC379" class="blob-code js-file-line"><span class="pl-c">// the number of dimensions to describe the location of an echo </span></td>
      </tr>
      <tr>
        <td id="L380" class="blob-num js-line-number" data-line-number="380"></td>
        <td id="LC380" class="blob-code js-file-line"><span class="pl-c">// (azimuth and elevAngle) as used in the &#39;pointsSelection&#39; array </span></td>
      </tr>
      <tr>
        <td id="L381" class="blob-num js-line-number" data-line-number="381"></td>
        <td id="LC381" class="blob-code js-file-line"><span class="pl-c">// that is passed to svdfit</span></td>
      </tr>
      <tr>
        <td id="L382" class="blob-num js-line-number" data-line-number="382"></td>
        <td id="LC382" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nDims;</td>
      </tr>
      <tr>
        <td id="L383" class="blob-num js-line-number" data-line-number="383"></td>
        <td id="LC383" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L384" class="blob-num js-line-number" data-line-number="384"></td>
        <td id="LC384" class="blob-code js-file-line"><span class="pl-c">// how many parameters are fitted by the svdfit procedure</span></td>
      </tr>
      <tr>
        <td id="L385" class="blob-num js-line-number" data-line-number="385"></td>
        <td id="LC385" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> nParsFitted;</td>
      </tr>
      <tr>
        <td id="L386" class="blob-num js-line-number" data-line-number="386"></td>
        <td id="LC386" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L387" class="blob-num js-line-number" data-line-number="387"></td>
        <td id="LC387" class="blob-code js-file-line"><span class="pl-c">// the factor that is used when converting from Z to eta</span></td>
      </tr>
      <tr>
        <td id="L388" class="blob-num js-line-number" data-line-number="388"></td>
        <td id="LC388" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> dbzFactor;</td>
      </tr>
      <tr>
        <td id="L389" class="blob-num js-line-number" data-line-number="389"></td>
        <td id="LC389" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L390" class="blob-num js-line-number" data-line-number="390"></td>
        <td id="LC390" class="blob-code js-file-line"><span class="pl-c">// whether the vol2bird module has been initialized</span></td>
      </tr>
      <tr>
        <td id="L391" class="blob-num js-line-number" data-line-number="391"></td>
        <td id="LC391" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> initializationSuccessful = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L392" class="blob-num js-line-number" data-line-number="392"></td>
        <td id="LC392" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L393" class="blob-num js-line-number" data-line-number="393"></td>
        <td id="LC393" class="blob-code js-file-line"><span class="pl-c">// during calculation of iProfileType == 3, use this array to store the</span></td>
      </tr>
      <tr>
        <td id="L394" class="blob-num js-line-number" data-line-number="394"></td>
        <td id="LC394" class="blob-code js-file-line"><span class="pl-c">// result of (chi &lt; stdDevMinBird), such that it can be used later during </span></td>
      </tr>
      <tr>
        <td id="L395" class="blob-num js-line-number" data-line-number="395"></td>
        <td id="LC395" class="blob-code js-file-line"><span class="pl-c">// calculation of iProfileType == 1</span></td>
      </tr>
      <tr>
        <td id="L396" class="blob-num js-line-number" data-line-number="396"></td>
        <td id="LC396" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span>* scatterersAreNotBirds;</td>
      </tr>
      <tr>
        <td id="L397" class="blob-num js-line-number" data-line-number="397"></td>
        <td id="LC397" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L398" class="blob-num js-line-number" data-line-number="398"></td>
        <td id="LC398" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L399" class="blob-num js-line-number" data-line-number="399"></td>
        <td id="LC399" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L400" class="blob-num js-line-number" data-line-number="400"></td>
        <td id="LC400" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L401" class="blob-num js-line-number" data-line-number="401"></td>
        <td id="LC401" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L402" class="blob-num js-line-number" data-line-number="402"></td>
        <td id="LC402" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">analyzeCells</span>(<span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *dbzImage, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *vradImage,</td>
      </tr>
      <tr>
        <td id="L403" class="blob-num js-line-number" data-line-number="403"></td>
        <td id="LC403" class="blob-code js-file-line">        <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *texImage, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *clutterImage, <span class="pl-k">int</span> *cellImage,</td>
      </tr>
      <tr>
        <td id="L404" class="blob-num js-line-number" data-line-number="404"></td>
        <td id="LC404" class="blob-code js-file-line">        <span class="pl-k">const</span> SCANMETA *dbzMeta, <span class="pl-k">const</span> SCANMETA *vradMeta, <span class="pl-k">const</span> SCANMETA *texMeta, <span class="pl-k">const</span> SCANMETA *clutterMeta,</td>
      </tr>
      <tr>
        <td id="L405" class="blob-num js-line-number" data-line-number="405"></td>
        <td id="LC405" class="blob-code js-file-line">        <span class="pl-k">const</span> <span class="pl-k">int</span> nCells, <span class="pl-k">const</span> <span class="pl-k">int</span> useStaticClutterData) {</td>
      </tr>
      <tr>
        <td id="L406" class="blob-num js-line-number" data-line-number="406"></td>
        <td id="LC406" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L407" class="blob-num js-line-number" data-line-number="407"></td>
        <td id="LC407" class="blob-code js-file-line">    <span class="pl-c">// ----------------------------------------------------------------------------------- // </span></td>
      </tr>
      <tr>
        <td id="L408" class="blob-num js-line-number" data-line-number="408"></td>
        <td id="LC408" class="blob-code js-file-line">    <span class="pl-c">//  This function analyzes the cellImage array found by the &#39;findWeatherCells&#39;         //</span></td>
      </tr>
      <tr>
        <td id="L409" class="blob-num js-line-number" data-line-number="409"></td>
        <td id="LC409" class="blob-code js-file-line">    <span class="pl-c">//  procedure. Small cells are rejected and the cells are re-numbered according        //</span></td>
      </tr>
      <tr>
        <td id="L410" class="blob-num js-line-number" data-line-number="410"></td>
        <td id="LC410" class="blob-code js-file-line">    <span class="pl-c">//  to size. The final number of cells in cellImage is returned as an integer.         // </span></td>
      </tr>
      <tr>
        <td id="L411" class="blob-num js-line-number" data-line-number="411"></td>
        <td id="LC411" class="blob-code js-file-line">    <span class="pl-c">// ----------------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L412" class="blob-num js-line-number" data-line-number="412"></td>
        <td id="LC412" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L413" class="blob-num js-line-number" data-line-number="413"></td>
        <td id="LC413" class="blob-code js-file-line">    CELLPROP *cellProp;</td>
      </tr>
      <tr>
        <td id="L414" class="blob-num js-line-number" data-line-number="414"></td>
        <td id="LC414" class="blob-code js-file-line">    <span class="pl-k">int</span> iCell;</td>
      </tr>
      <tr>
        <td id="L415" class="blob-num js-line-number" data-line-number="415"></td>
        <td id="LC415" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L416" class="blob-num js-line-number" data-line-number="416"></td>
        <td id="LC416" class="blob-code js-file-line">    <span class="pl-k">int</span> nGlobal;</td>
      </tr>
      <tr>
        <td id="L417" class="blob-num js-line-number" data-line-number="417"></td>
        <td id="LC417" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L418" class="blob-num js-line-number" data-line-number="418"></td>
        <td id="LC418" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzim;</td>
      </tr>
      <tr>
        <td id="L419" class="blob-num js-line-number" data-line-number="419"></td>
        <td id="LC419" class="blob-code js-file-line">    <span class="pl-k">int</span> nCellsValid;</td>
      </tr>
      <tr>
        <td id="L420" class="blob-num js-line-number" data-line-number="420"></td>
        <td id="LC420" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzim;</td>
      </tr>
      <tr>
        <td id="L421" class="blob-num js-line-number" data-line-number="421"></td>
        <td id="LC421" class="blob-code js-file-line">    <span class="pl-k">int</span> nRang;</td>
      </tr>
      <tr>
        <td id="L422" class="blob-num js-line-number" data-line-number="422"></td>
        <td id="LC422" class="blob-code js-file-line">    <span class="pl-k">float</span> validArea;</td>
      </tr>
      <tr>
        <td id="L423" class="blob-num js-line-number" data-line-number="423"></td>
        <td id="LC423" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzValue;</td>
      </tr>
      <tr>
        <td id="L424" class="blob-num js-line-number" data-line-number="424"></td>
        <td id="LC424" class="blob-code js-file-line">    <span class="pl-k">float</span> texValue;</td>
      </tr>
      <tr>
        <td id="L425" class="blob-num js-line-number" data-line-number="425"></td>
        <td id="LC425" class="blob-code js-file-line">    <span class="pl-k">float</span> vradValue;</td>
      </tr>
      <tr>
        <td id="L426" class="blob-num js-line-number" data-line-number="426"></td>
        <td id="LC426" class="blob-code js-file-line">    <span class="pl-k">float</span> clutterValue;</td>
      </tr>
      <tr>
        <td id="L427" class="blob-num js-line-number" data-line-number="427"></td>
        <td id="LC427" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L428" class="blob-num js-line-number" data-line-number="428"></td>
        <td id="LC428" class="blob-code js-file-line">    nCellsValid = nCells;</td>
      </tr>
      <tr>
        <td id="L429" class="blob-num js-line-number" data-line-number="429"></td>
        <td id="LC429" class="blob-code js-file-line">    nRang = dbzMeta-&gt;nRang;</td>
      </tr>
      <tr>
        <td id="L430" class="blob-num js-line-number" data-line-number="430"></td>
        <td id="LC430" class="blob-code js-file-line">    nAzim = dbzMeta-&gt;nAzim;</td>
      </tr>
      <tr>
        <td id="L431" class="blob-num js-line-number" data-line-number="431"></td>
        <td id="LC431" class="blob-code js-file-line">    nGlobal = nAzim*nRang;</td>
      </tr>
      <tr>
        <td id="L432" class="blob-num js-line-number" data-line-number="432"></td>
        <td id="LC432" class="blob-code js-file-line">    nCellsValid = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L433" class="blob-num js-line-number" data-line-number="433"></td>
        <td id="LC433" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L434" class="blob-num js-line-number" data-line-number="434"></td>
        <td id="LC434" class="blob-code js-file-line">    <span class="pl-k">if</span> (nCells == <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L435" class="blob-num js-line-number" data-line-number="435"></td>
        <td id="LC435" class="blob-code js-file-line">        iGlobal = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L436" class="blob-num js-line-number" data-line-number="436"></td>
        <td id="LC436" class="blob-code js-file-line">        <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) {</td>
      </tr>
      <tr>
        <td id="L437" class="blob-num js-line-number" data-line-number="437"></td>
        <td id="LC437" class="blob-code js-file-line">            <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L438" class="blob-num js-line-number" data-line-number="438"></td>
        <td id="LC438" class="blob-code js-file-line">                cellImage[iGlobal] = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L439" class="blob-num js-line-number" data-line-number="439"></td>
        <td id="LC439" class="blob-code js-file-line">                iGlobal++;</td>
      </tr>
      <tr>
        <td id="L440" class="blob-num js-line-number" data-line-number="440"></td>
        <td id="LC440" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L441" class="blob-num js-line-number" data-line-number="441"></td>
        <td id="LC441" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L442" class="blob-num js-line-number" data-line-number="442"></td>
        <td id="LC442" class="blob-code js-file-line">        <span class="pl-k">return</span> nCellsValid;</td>
      </tr>
      <tr>
        <td id="L443" class="blob-num js-line-number" data-line-number="443"></td>
        <td id="LC443" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L444" class="blob-num js-line-number" data-line-number="444"></td>
        <td id="LC444" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L445" class="blob-num js-line-number" data-line-number="445"></td>
        <td id="LC445" class="blob-code js-file-line">    <span class="pl-c">// Allocating and initializing memory for cell properties.</span></td>
      </tr>
      <tr>
        <td id="L446" class="blob-num js-line-number" data-line-number="446"></td>
        <td id="LC446" class="blob-code js-file-line">    cellProp = (CELLPROP *)<span class="pl-c1">malloc</span>(nCells*<span class="pl-k">sizeof</span>(CELLPROP));</td>
      </tr>
      <tr>
        <td id="L447" class="blob-num js-line-number" data-line-number="447"></td>
        <td id="LC447" class="blob-code js-file-line">    <span class="pl-k">if</span> (!cellProp) {</td>
      </tr>
      <tr>
        <td id="L448" class="blob-num js-line-number" data-line-number="448"></td>
        <td id="LC448" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Requested memory could not be allocated!<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L449" class="blob-num js-line-number" data-line-number="449"></td>
        <td id="LC449" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">10</span>;</td>
      </tr>
      <tr>
        <td id="L450" class="blob-num js-line-number" data-line-number="450"></td>
        <td id="LC450" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L451" class="blob-num js-line-number" data-line-number="451"></td>
        <td id="LC451" class="blob-code js-file-line">    <span class="pl-k">for</span> (iCell = <span class="pl-c1">0</span>; iCell &lt; nCells; iCell++) {</td>
      </tr>
      <tr>
        <td id="L452" class="blob-num js-line-number" data-line-number="452"></td>
        <td id="LC452" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">iRangOfMax</span> = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L453" class="blob-num js-line-number" data-line-number="453"></td>
        <td id="LC453" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">iAzimOfMax</span> = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L454" class="blob-num js-line-number" data-line-number="454"></td>
        <td id="LC454" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">nGates</span> = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L455" class="blob-num js-line-number" data-line-number="455"></td>
        <td id="LC455" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">nGatesClutter</span> = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L456" class="blob-num js-line-number" data-line-number="456"></td>
        <td id="LC456" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">dbzAvg</span> = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L457" class="blob-num js-line-number" data-line-number="457"></td>
        <td id="LC457" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">texAvg</span> = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L458" class="blob-num js-line-number" data-line-number="458"></td>
        <td id="LC458" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">dbzMax</span> = dbzMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L459" class="blob-num js-line-number" data-line-number="459"></td>
        <td id="LC459" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">index</span> = iCell;</td>
      </tr>
      <tr>
        <td id="L460" class="blob-num js-line-number" data-line-number="460"></td>
        <td id="LC460" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">drop</span> = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L461" class="blob-num js-line-number" data-line-number="461"></td>
        <td id="LC461" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">cv</span> = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L462" class="blob-num js-line-number" data-line-number="462"></td>
        <td id="LC462" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L463" class="blob-num js-line-number" data-line-number="463"></td>
        <td id="LC463" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L464" class="blob-num js-line-number" data-line-number="464"></td>
        <td id="LC464" class="blob-code js-file-line">    <span class="pl-c">// Calculation of cell properties.</span></td>
      </tr>
      <tr>
        <td id="L465" class="blob-num js-line-number" data-line-number="465"></td>
        <td id="LC465" class="blob-code js-file-line">    <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) {</td>
      </tr>
      <tr>
        <td id="L466" class="blob-num js-line-number" data-line-number="466"></td>
        <td id="LC466" class="blob-code js-file-line">        <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L467" class="blob-num js-line-number" data-line-number="467"></td>
        <td id="LC467" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L468" class="blob-num js-line-number" data-line-number="468"></td>
        <td id="LC468" class="blob-code js-file-line">            iGlobal = iRang + iAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L469" class="blob-num js-line-number" data-line-number="469"></td>
        <td id="LC469" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L470" class="blob-num js-line-number" data-line-number="470"></td>
        <td id="LC470" class="blob-code js-file-line">            dbzValue = dbzMeta-&gt;valueScale * (<span class="pl-k">float</span>) dbzImage[iGlobal] + dbzMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L471" class="blob-num js-line-number" data-line-number="471"></td>
        <td id="LC471" class="blob-code js-file-line">            vradValue = vradMeta-&gt;valueScale * (<span class="pl-k">float</span>) vradImage[iGlobal] + vradMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L472" class="blob-num js-line-number" data-line-number="472"></td>
        <td id="LC472" class="blob-code js-file-line">            clutterValue = clutterMeta-&gt;valueScale * (<span class="pl-k">float</span>) clutterImage[iGlobal] + clutterMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L473" class="blob-num js-line-number" data-line-number="473"></td>
        <td id="LC473" class="blob-code js-file-line">            texValue = texMeta-&gt;valueScale * (<span class="pl-k">float</span>) texImage[iGlobal] + texMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L474" class="blob-num js-line-number" data-line-number="474"></td>
        <td id="LC474" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L475" class="blob-num js-line-number" data-line-number="475"></td>
        <td id="LC475" class="blob-code js-file-line">            iCell = cellImage[iGlobal];</td>
      </tr>
      <tr>
        <td id="L476" class="blob-num js-line-number" data-line-number="476"></td>
        <td id="LC476" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L477" class="blob-num js-line-number" data-line-number="477"></td>
        <td id="LC477" class="blob-code js-file-line">            <span class="pl-k">if</span> (iCell&lt;<span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L478" class="blob-num js-line-number" data-line-number="478"></td>
        <td id="LC478" class="blob-code js-file-line">                <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L479" class="blob-num js-line-number" data-line-number="479"></td>
        <td id="LC479" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L480" class="blob-num js-line-number" data-line-number="480"></td>
        <td id="LC480" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L481" class="blob-num js-line-number" data-line-number="481"></td>
        <td id="LC481" class="blob-code js-file-line">            #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L482" class="blob-num js-line-number" data-line-number="482"></td>
        <td id="LC482" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>dbzValue = <span class="pl-c1">%f</span>; vradValue = <span class="pl-c1">%f</span>; clutterValue = <span class="pl-c1">%f</span>; texValue = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,dbzValue,vradValue,clutterValue,texValue);</td>
      </tr>
      <tr>
        <td id="L483" class="blob-num js-line-number" data-line-number="483"></td>
        <td id="LC483" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>iGlobal = <span class="pl-c1">%d</span>, iCell = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iGlobal,iCell);</td>
      </tr>
      <tr>
        <td id="L484" class="blob-num js-line-number" data-line-number="484"></td>
        <td id="LC484" class="blob-code js-file-line">            #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L485" class="blob-num js-line-number" data-line-number="485"></td>
        <td id="LC485" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L486" class="blob-num js-line-number" data-line-number="486"></td>
        <td id="LC486" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">nGates</span> += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L487" class="blob-num js-line-number" data-line-number="487"></td>
        <td id="LC487" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L488" class="blob-num js-line-number" data-line-number="488"></td>
        <td id="LC488" class="blob-code js-file-line">            <span class="pl-c">// low radial velocities are treated as clutter, not included in calculation cell properties</span></td>
      </tr>
      <tr>
        <td id="L489" class="blob-num js-line-number" data-line-number="489"></td>
        <td id="LC489" class="blob-code js-file-line">            <span class="pl-k">if</span> (vradValue &lt; vradMin){</td>
      </tr>
      <tr>
        <td id="L490" class="blob-num js-line-number" data-line-number="490"></td>
        <td id="LC490" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L491" class="blob-num js-line-number" data-line-number="491"></td>
        <td id="LC491" class="blob-code js-file-line">                cellProp[iCell].<span class="pl-smi">nGatesClutter</span> += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L492" class="blob-num js-line-number" data-line-number="492"></td>
        <td id="LC492" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L493" class="blob-num js-line-number" data-line-number="493"></td>
        <td id="LC493" class="blob-code js-file-line">                #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L494" class="blob-num js-line-number" data-line-number="494"></td>
        <td id="LC494" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>iGlobal = <span class="pl-c1">%d</span>: vrad too low...treating as clutter<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iGlobal);</td>
      </tr>
      <tr>
        <td id="L495" class="blob-num js-line-number" data-line-number="495"></td>
        <td id="LC495" class="blob-code js-file-line">                #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L496" class="blob-num js-line-number" data-line-number="496"></td>
        <td id="LC496" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L497" class="blob-num js-line-number" data-line-number="497"></td>
        <td id="LC497" class="blob-code js-file-line">                <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L498" class="blob-num js-line-number" data-line-number="498"></td>
        <td id="LC498" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L499" class="blob-num js-line-number" data-line-number="499"></td>
        <td id="LC499" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L500" class="blob-num js-line-number" data-line-number="500"></td>
        <td id="LC500" class="blob-code js-file-line">            <span class="pl-c">// pixels in clutter map not included in calculation cell properties</span></td>
      </tr>
      <tr>
        <td id="L501" class="blob-num js-line-number" data-line-number="501"></td>
        <td id="LC501" class="blob-code js-file-line">            <span class="pl-k">if</span> (useStaticClutterData == <span class="pl-c1">TRUE</span>){</td>
      </tr>
      <tr>
        <td id="L502" class="blob-num js-line-number" data-line-number="502"></td>
        <td id="LC502" class="blob-code js-file-line">                <span class="pl-k">if</span> (clutterValue &gt; clutterValueMin){</td>
      </tr>
      <tr>
        <td id="L503" class="blob-num js-line-number" data-line-number="503"></td>
        <td id="LC503" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">nGatesClutter</span> += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L504" class="blob-num js-line-number" data-line-number="504"></td>
        <td id="LC504" class="blob-code js-file-line">                    <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L505" class="blob-num js-line-number" data-line-number="505"></td>
        <td id="LC505" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L506" class="blob-num js-line-number" data-line-number="506"></td>
        <td id="LC506" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L507" class="blob-num js-line-number" data-line-number="507"></td>
        <td id="LC507" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L508" class="blob-num js-line-number" data-line-number="508"></td>
        <td id="LC508" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L509" class="blob-num js-line-number" data-line-number="509"></td>
        <td id="LC509" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L510" class="blob-num js-line-number" data-line-number="510"></td>
        <td id="LC510" class="blob-code js-file-line">            <span class="pl-k">if</span> (dbzValue &gt; cellProp[iCell].<span class="pl-smi">dbzMax</span>) {</td>
      </tr>
      <tr>
        <td id="L511" class="blob-num js-line-number" data-line-number="511"></td>
        <td id="LC511" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L512" class="blob-num js-line-number" data-line-number="512"></td>
        <td id="LC512" class="blob-code js-file-line">                #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L513" class="blob-num js-line-number" data-line-number="513"></td>
        <td id="LC513" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%d</span>: new dbzMax value of <span class="pl-c1">%f</span> found for this cell (<span class="pl-c1">%d</span>).<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iGlobal,dbzValue,iCell);</td>
      </tr>
      <tr>
        <td id="L514" class="blob-num js-line-number" data-line-number="514"></td>
        <td id="LC514" class="blob-code js-file-line">                #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L515" class="blob-num js-line-number" data-line-number="515"></td>
        <td id="LC515" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L516" class="blob-num js-line-number" data-line-number="516"></td>
        <td id="LC516" class="blob-code js-file-line">                cellProp[iCell].<span class="pl-smi">dbzMax</span> = dbzValue;</td>
      </tr>
      <tr>
        <td id="L517" class="blob-num js-line-number" data-line-number="517"></td>
        <td id="LC517" class="blob-code js-file-line">                cellProp[iCell].<span class="pl-smi">iRangOfMax</span> = iGlobal%nRang;</td>
      </tr>
      <tr>
        <td id="L518" class="blob-num js-line-number" data-line-number="518"></td>
        <td id="LC518" class="blob-code js-file-line">                cellProp[iCell].<span class="pl-smi">iAzimOfMax</span> = iGlobal/nRang;</td>
      </tr>
      <tr>
        <td id="L519" class="blob-num js-line-number" data-line-number="519"></td>
        <td id="LC519" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L520" class="blob-num js-line-number" data-line-number="520"></td>
        <td id="LC520" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">dbzAvg</span> += dbzValue;</td>
      </tr>
      <tr>
        <td id="L521" class="blob-num js-line-number" data-line-number="521"></td>
        <td id="LC521" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">texAvg</span> += texValue;</td>
      </tr>
      <tr>
        <td id="L522" class="blob-num js-line-number" data-line-number="522"></td>
        <td id="LC522" class="blob-code js-file-line">        } <span class="pl-c">// for (iRang = 0; iRang &lt; nRang; iRang++)</span></td>
      </tr>
      <tr>
        <td id="L523" class="blob-num js-line-number" data-line-number="523"></td>
        <td id="LC523" class="blob-code js-file-line">    } <span class="pl-c">// for (iAzim = 0; iAzim &lt; nAzim; iAzim++)</span></td>
      </tr>
      <tr>
        <td id="L524" class="blob-num js-line-number" data-line-number="524"></td>
        <td id="LC524" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L525" class="blob-num js-line-number" data-line-number="525"></td>
        <td id="LC525" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L526" class="blob-num js-line-number" data-line-number="526"></td>
        <td id="LC526" class="blob-code js-file-line">    <span class="pl-k">for</span> (iCell = <span class="pl-c1">0</span>; iCell &lt; nCells; iCell++) {</td>
      </tr>
      <tr>
        <td id="L527" class="blob-num js-line-number" data-line-number="527"></td>
        <td id="LC527" class="blob-code js-file-line">        validArea = cellProp[iCell].<span class="pl-smi">nGates</span> - cellProp[iCell].<span class="pl-smi">nGatesClutter</span>;</td>
      </tr>
      <tr>
        <td id="L528" class="blob-num js-line-number" data-line-number="528"></td>
        <td id="LC528" class="blob-code js-file-line">        <span class="pl-k">if</span> (validArea &gt; <span class="pl-c1">0</span>){</td>
      </tr>
      <tr>
        <td id="L529" class="blob-num js-line-number" data-line-number="529"></td>
        <td id="LC529" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">dbzAvg</span> /= validArea;</td>
      </tr>
      <tr>
        <td id="L530" class="blob-num js-line-number" data-line-number="530"></td>
        <td id="LC530" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">texAvg</span> /= validArea;</td>
      </tr>
      <tr>
        <td id="L531" class="blob-num js-line-number" data-line-number="531"></td>
        <td id="LC531" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">cv</span> = cellProp[iCell].<span class="pl-smi">texAvg</span> / cellProp[iCell].<span class="pl-smi">dbzAvg</span>;</td>
      </tr>
      <tr>
        <td id="L532" class="blob-num js-line-number" data-line-number="532"></td>
        <td id="LC532" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L533" class="blob-num js-line-number" data-line-number="533"></td>
        <td id="LC533" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L534" class="blob-num js-line-number" data-line-number="534"></td>
        <td id="LC534" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L535" class="blob-num js-line-number" data-line-number="535"></td>
        <td id="LC535" class="blob-code js-file-line">    <span class="pl-c">// determine which cells to drop from map based on low mean dBZ / high stdev /</span></td>
      </tr>
      <tr>
        <td id="L536" class="blob-num js-line-number" data-line-number="536"></td>
        <td id="LC536" class="blob-code js-file-line">    <span class="pl-c">// small area / high percentage clutter</span></td>
      </tr>
      <tr>
        <td id="L537" class="blob-num js-line-number" data-line-number="537"></td>
        <td id="LC537" class="blob-code js-file-line">    <span class="pl-k">for</span> (iCell = <span class="pl-c1">0</span>; iCell &lt; nCells; iCell++) {</td>
      </tr>
      <tr>
        <td id="L538" class="blob-num js-line-number" data-line-number="538"></td>
        <td id="LC538" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellProp[iCell].<span class="pl-smi">nGates</span> &lt; nGatesCellMin ||</td>
      </tr>
      <tr>
        <td id="L539" class="blob-num js-line-number" data-line-number="539"></td>
        <td id="LC539" class="blob-code js-file-line">            (cellProp[iCell].<span class="pl-smi">dbzAvg</span> &lt; cellDbzMin &amp;&amp;</td>
      </tr>
      <tr>
        <td id="L540" class="blob-num js-line-number" data-line-number="540"></td>
        <td id="LC540" class="blob-code js-file-line">             cellProp[iCell].<span class="pl-smi">texAvg</span> &gt; cellStdDevMax &amp;&amp;</td>
      </tr>
      <tr>
        <td id="L541" class="blob-num js-line-number" data-line-number="541"></td>
        <td id="LC541" class="blob-code js-file-line">             ((<span class="pl-k">float</span>) cellProp[iCell].<span class="pl-smi">nGatesClutter</span> / cellProp[iCell].<span class="pl-smi">nGates</span>) &lt; cellClutterFractionMax )) { <span class="pl-c">// FIXME this condition still looks wrong if you ask me</span></td>
      </tr>
      <tr>
        <td id="L542" class="blob-num js-line-number" data-line-number="542"></td>
        <td id="LC542" class="blob-code js-file-line">            <span class="pl-c">// Terms 2,3 and 4 are combined with &amp;&amp; to be conservative in labeling stuff as</span></td>
      </tr>
      <tr>
        <td id="L543" class="blob-num js-line-number" data-line-number="543"></td>
        <td id="LC543" class="blob-code js-file-line">            <span class="pl-c">// bird migration --see discussion of issue #37 on GitHub.</span></td>
      </tr>
      <tr>
        <td id="L544" class="blob-num js-line-number" data-line-number="544"></td>
        <td id="LC544" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">drop</span> = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L545" class="blob-num js-line-number" data-line-number="545"></td>
        <td id="LC545" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">nGates</span> = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L546" class="blob-num js-line-number" data-line-number="546"></td>
        <td id="LC546" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L547" class="blob-num js-line-number" data-line-number="547"></td>
        <td id="LC547" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L548" class="blob-num js-line-number" data-line-number="548"></td>
        <td id="LC548" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L549" class="blob-num js-line-number" data-line-number="549"></td>
        <td id="LC549" class="blob-code js-file-line">    <span class="pl-c">// sorting cell properties according to cell area. Drop small cells from map</span></td>
      </tr>
      <tr>
        <td id="L550" class="blob-num js-line-number" data-line-number="550"></td>
        <td id="LC550" class="blob-code js-file-line">    nCellsValid = <span class="pl-c1">updateMap</span>(cellImage,nGlobal,cellProp,nCells);</td>
      </tr>
      <tr>
        <td id="L551" class="blob-num js-line-number" data-line-number="551"></td>
        <td id="LC551" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L552" class="blob-num js-line-number" data-line-number="552"></td>
        <td id="LC552" class="blob-code js-file-line">    <span class="pl-c">// printing of cell properties to stderr</span></td>
      </tr>
      <tr>
        <td id="L553" class="blob-num js-line-number" data-line-number="553"></td>
        <td id="LC553" class="blob-code js-file-line">    <span class="pl-k">if</span> (printCellProp == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L554" class="blob-num js-line-number" data-line-number="554"></td>
        <td id="LC554" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>#Cell analysis for elevation <span class="pl-c1">%f</span>:<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,dbzMeta-&gt;elev);</td>
      </tr>
      <tr>
        <td id="L555" class="blob-num js-line-number" data-line-number="555"></td>
        <td id="LC555" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>#Minimum cell area in pixels   : <span class="pl-c1">%i</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,nGatesCellMin);</td>
      </tr>
      <tr>
        <td id="L556" class="blob-num js-line-number" data-line-number="556"></td>
        <td id="LC556" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>#Threshold for mean dBZ cell   : <span class="pl-c1">%g</span> dBZ<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,cellDbzMin);</td>
      </tr>
      <tr>
        <td id="L557" class="blob-num js-line-number" data-line-number="557"></td>
        <td id="LC557" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>#Threshold for mean stdev cell : <span class="pl-c1">%g</span> dBZ<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,cellStdDevMax);</td>
      </tr>
      <tr>
        <td id="L558" class="blob-num js-line-number" data-line-number="558"></td>
        <td id="LC558" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>#Valid cells                   : <span class="pl-c1">%i</span>/<span class="pl-c1">%i</span><span class="pl-cce">\n</span>#<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,nCellsValid,nCells);</td>
      </tr>
      <tr>
        <td id="L559" class="blob-num js-line-number" data-line-number="559"></td>
        <td id="LC559" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>cellProp: .index .nGates .nGatesClutter .dbzAvg .texAvg .cv   .dbzMax .iRangOfMax .iAzimOfMax .drop<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L560" class="blob-num js-line-number" data-line-number="560"></td>
        <td id="LC560" class="blob-code js-file-line">        <span class="pl-k">for</span> (iCell = <span class="pl-c1">0</span>; iCell &lt; nCells; iCell++) {</td>
      </tr>
      <tr>
        <td id="L561" class="blob-num js-line-number" data-line-number="561"></td>
        <td id="LC561" class="blob-code js-file-line">            <span class="pl-k">if</span> (cellProp[iCell].<span class="pl-smi">nGates</span> == <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L562" class="blob-num js-line-number" data-line-number="562"></td>
        <td id="LC562" class="blob-code js-file-line">                <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L563" class="blob-num js-line-number" data-line-number="563"></td>
        <td id="LC563" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L564" class="blob-num js-line-number" data-line-number="564"></td>
        <td id="LC564" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>cellProp: <span class="pl-c1">%6d</span> <span class="pl-c1">%7d</span> <span class="pl-c1">%14d</span> <span class="pl-c1">%7.2f</span> <span class="pl-c1">%7.2f</span> <span class="pl-c1">%5.2f</span> <span class="pl-c1">%7.2f</span> <span class="pl-c1">%11d</span> <span class="pl-c1">%11d</span> <span class="pl-c1">%5c</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L565" class="blob-num js-line-number" data-line-number="565"></td>
        <td id="LC565" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">index</span>,</td>
      </tr>
      <tr>
        <td id="L566" class="blob-num js-line-number" data-line-number="566"></td>
        <td id="LC566" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">nGates</span>,</td>
      </tr>
      <tr>
        <td id="L567" class="blob-num js-line-number" data-line-number="567"></td>
        <td id="LC567" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">nGatesClutter</span>,</td>
      </tr>
      <tr>
        <td id="L568" class="blob-num js-line-number" data-line-number="568"></td>
        <td id="LC568" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">dbzAvg</span>,</td>
      </tr>
      <tr>
        <td id="L569" class="blob-num js-line-number" data-line-number="569"></td>
        <td id="LC569" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">texAvg</span>,</td>
      </tr>
      <tr>
        <td id="L570" class="blob-num js-line-number" data-line-number="570"></td>
        <td id="LC570" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">cv</span>,</td>
      </tr>
      <tr>
        <td id="L571" class="blob-num js-line-number" data-line-number="571"></td>
        <td id="LC571" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">dbzMax</span>,</td>
      </tr>
      <tr>
        <td id="L572" class="blob-num js-line-number" data-line-number="572"></td>
        <td id="LC572" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">iRangOfMax</span>,</td>
      </tr>
      <tr>
        <td id="L573" class="blob-num js-line-number" data-line-number="573"></td>
        <td id="LC573" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">iAzimOfMax</span>,</td>
      </tr>
      <tr>
        <td id="L574" class="blob-num js-line-number" data-line-number="574"></td>
        <td id="LC574" class="blob-code js-file-line">                    cellProp[iCell].<span class="pl-smi">drop</span> == <span class="pl-c1">TRUE</span> ? <span class="pl-s"><span class="pl-pds">&#39;</span>T<span class="pl-pds">&#39;</span></span> : <span class="pl-s"><span class="pl-pds">&#39;</span>F<span class="pl-pds">&#39;</span></span>);</td>
      </tr>
      <tr>
        <td id="L575" class="blob-num js-line-number" data-line-number="575"></td>
        <td id="LC575" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L576" class="blob-num js-line-number" data-line-number="576"></td>
        <td id="LC576" class="blob-code js-file-line">    } <span class="pl-c">// endif (printCellProp == TRUE)</span></td>
      </tr>
      <tr>
        <td id="L577" class="blob-num js-line-number" data-line-number="577"></td>
        <td id="LC577" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L578" class="blob-num js-line-number" data-line-number="578"></td>
        <td id="LC578" class="blob-code js-file-line">    <span class="pl-c1">free</span>(cellProp);</td>
      </tr>
      <tr>
        <td id="L579" class="blob-num js-line-number" data-line-number="579"></td>
        <td id="LC579" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L580" class="blob-num js-line-number" data-line-number="580"></td>
        <td id="LC580" class="blob-code js-file-line">    <span class="pl-k">return</span> nCellsValid;</td>
      </tr>
      <tr>
        <td id="L581" class="blob-num js-line-number" data-line-number="581"></td>
        <td id="LC581" class="blob-code js-file-line">} <span class="pl-c">// analyzeCells</span></td>
      </tr>
      <tr>
        <td id="L582" class="blob-num js-line-number" data-line-number="582"></td>
        <td id="LC582" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L583" class="blob-num js-line-number" data-line-number="583"></td>
        <td id="LC583" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L584" class="blob-num js-line-number" data-line-number="584"></td>
        <td id="LC584" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L585" class="blob-num js-line-number" data-line-number="585"></td>
        <td id="LC585" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L586" class="blob-num js-line-number" data-line-number="586"></td>
        <td id="LC586" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L587" class="blob-num js-line-number" data-line-number="587"></td>
        <td id="LC587" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">float</span> <span class="pl-en">calcDist</span>(<span class="pl-k">const</span> <span class="pl-k">int</span> iRang1, <span class="pl-k">const</span> <span class="pl-k">int</span> iAzim1, <span class="pl-k">const</span> <span class="pl-k">int</span> iRang2, <span class="pl-k">const</span> </td>
      </tr>
      <tr>
        <td id="L588" class="blob-num js-line-number" data-line-number="588"></td>
        <td id="LC588" class="blob-code js-file-line">    <span class="pl-k">const</span> <span class="pl-k">int</span> iAzim2, <span class="pl-k">const</span> <span class="pl-k">float</span> rangScale, <span class="pl-k">const</span> <span class="pl-k">float</span> azimScaleDeg) {</td>
      </tr>
      <tr>
        <td id="L589" class="blob-num js-line-number" data-line-number="589"></td>
        <td id="LC589" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L590" class="blob-num js-line-number" data-line-number="590"></td>
        <td id="LC590" class="blob-code js-file-line">    <span class="pl-c">// -------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L591" class="blob-num js-line-number" data-line-number="591"></td>
        <td id="LC591" class="blob-code js-file-line">    <span class="pl-c">// This function calculates the distance between two gates  //</span></td>
      </tr>
      <tr>
        <td id="L592" class="blob-num js-line-number" data-line-number="592"></td>
        <td id="LC592" class="blob-code js-file-line">    <span class="pl-c">// -------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L593" class="blob-num js-line-number" data-line-number="593"></td>
        <td id="LC593" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L594" class="blob-num js-line-number" data-line-number="594"></td>
        <td id="LC594" class="blob-code js-file-line">    <span class="pl-k">float</span> range1;</td>
      </tr>
      <tr>
        <td id="L595" class="blob-num js-line-number" data-line-number="595"></td>
        <td id="LC595" class="blob-code js-file-line">    <span class="pl-k">float</span> range2;</td>
      </tr>
      <tr>
        <td id="L596" class="blob-num js-line-number" data-line-number="596"></td>
        <td id="LC596" class="blob-code js-file-line">    <span class="pl-k">float</span> azimuth1;</td>
      </tr>
      <tr>
        <td id="L597" class="blob-num js-line-number" data-line-number="597"></td>
        <td id="LC597" class="blob-code js-file-line">    <span class="pl-k">float</span> azimuth2;</td>
      </tr>
      <tr>
        <td id="L598" class="blob-num js-line-number" data-line-number="598"></td>
        <td id="LC598" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L599" class="blob-num js-line-number" data-line-number="599"></td>
        <td id="LC599" class="blob-code js-file-line">    range1 = (<span class="pl-k">float</span>) iRang1 * rangScale;</td>
      </tr>
      <tr>
        <td id="L600" class="blob-num js-line-number" data-line-number="600"></td>
        <td id="LC600" class="blob-code js-file-line">    range2 = (<span class="pl-k">float</span>) iRang2 * rangScale;</td>
      </tr>
      <tr>
        <td id="L601" class="blob-num js-line-number" data-line-number="601"></td>
        <td id="LC601" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L602" class="blob-num js-line-number" data-line-number="602"></td>
        <td id="LC602" class="blob-code js-file-line">    azimuth1 = (<span class="pl-k">float</span>) iAzim1 * azimScaleDeg * (<span class="pl-k">float</span>) DEG2RAD;</td>
      </tr>
      <tr>
        <td id="L603" class="blob-num js-line-number" data-line-number="603"></td>
        <td id="LC603" class="blob-code js-file-line">    azimuth2 = (<span class="pl-k">float</span>) iAzim2 * azimScaleDeg * (<span class="pl-k">float</span>) DEG2RAD;</td>
      </tr>
      <tr>
        <td id="L604" class="blob-num js-line-number" data-line-number="604"></td>
        <td id="LC604" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L605" class="blob-num js-line-number" data-line-number="605"></td>
        <td id="LC605" class="blob-code js-file-line">    <span class="pl-k">return</span> (<span class="pl-k">float</span>) <span class="pl-c1">sqrt</span>((range1 * range1) +</td>
      </tr>
      <tr>
        <td id="L606" class="blob-num js-line-number" data-line-number="606"></td>
        <td id="LC606" class="blob-code js-file-line">                (range2 * range2) -</td>
      </tr>
      <tr>
        <td id="L607" class="blob-num js-line-number" data-line-number="607"></td>
        <td id="LC607" class="blob-code js-file-line">                <span class="pl-c1">2</span> * (range1 * range2) * <span class="pl-c1">cos</span>(azimuth1-azimuth2));</td>
      </tr>
      <tr>
        <td id="L608" class="blob-num js-line-number" data-line-number="608"></td>
        <td id="LC608" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L609" class="blob-num js-line-number" data-line-number="609"></td>
        <td id="LC609" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L610" class="blob-num js-line-number" data-line-number="610"></td>
        <td id="LC610" class="blob-code js-file-line">} <span class="pl-c">// calcDist</span></td>
      </tr>
      <tr>
        <td id="L611" class="blob-num js-line-number" data-line-number="611"></td>
        <td id="LC611" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L612" class="blob-num js-line-number" data-line-number="612"></td>
        <td id="LC612" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L613" class="blob-num js-line-number" data-line-number="613"></td>
        <td id="LC613" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L614" class="blob-num js-line-number" data-line-number="614"></td>
        <td id="LC614" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L615" class="blob-num js-line-number" data-line-number="615"></td>
        <td id="LC615" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L616" class="blob-num js-line-number" data-line-number="616"></td>
        <td id="LC616" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">calcTexture</span>(<span class="pl-k">unsigned</span> <span class="pl-k">char</span>* texImage, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* vradImage,</td>
      </tr>
      <tr>
        <td id="L617" class="blob-num js-line-number" data-line-number="617"></td>
        <td id="LC617" class="blob-code js-file-line">        <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* dbzImage, <span class="pl-k">const</span> SCANMETA* texMeta, <span class="pl-k">const</span> SCANMETA* vradMeta,</td>
      </tr>
      <tr>
        <td id="L618" class="blob-num js-line-number" data-line-number="618"></td>
        <td id="LC618" class="blob-code js-file-line">        <span class="pl-k">const</span> SCANMETA* dbzMeta) {</td>
      </tr>
      <tr>
        <td id="L619" class="blob-num js-line-number" data-line-number="619"></td>
        <td id="LC619" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L620" class="blob-num js-line-number" data-line-number="620"></td>
        <td id="LC620" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L621" class="blob-num js-line-number" data-line-number="621"></td>
        <td id="LC621" class="blob-code js-file-line">    <span class="pl-c">// --------------------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L622" class="blob-num js-line-number" data-line-number="622"></td>
        <td id="LC622" class="blob-code js-file-line">    <span class="pl-c">// This function computes a texture parameter based on a block of (nRangNeighborhood x     //</span></td>
      </tr>
      <tr>
        <td id="L623" class="blob-num js-line-number" data-line-number="623"></td>
        <td id="LC623" class="blob-code js-file-line">    <span class="pl-c">// nAzimNeighborhood) pixels. The texture parameter equals the local standard deviation    //</span></td>
      </tr>
      <tr>
        <td id="L624" class="blob-num js-line-number" data-line-number="624"></td>
        <td id="LC624" class="blob-code js-file-line">    <span class="pl-c">// in the radial velocity field                                                            //</span></td>
      </tr>
      <tr>
        <td id="L625" class="blob-num js-line-number" data-line-number="625"></td>
        <td id="LC625" class="blob-code js-file-line">    <span class="pl-c">// --------------------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L626" class="blob-num js-line-number" data-line-number="626"></td>
        <td id="LC626" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L627" class="blob-num js-line-number" data-line-number="627"></td>
        <td id="LC627" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L628" class="blob-num js-line-number" data-line-number="628"></td>
        <td id="LC628" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L629" class="blob-num js-line-number" data-line-number="629"></td>
        <td id="LC629" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzim;</td>
      </tr>
      <tr>
        <td id="L630" class="blob-num js-line-number" data-line-number="630"></td>
        <td id="LC630" class="blob-code js-file-line">    <span class="pl-k">int</span> nRang;</td>
      </tr>
      <tr>
        <td id="L631" class="blob-num js-line-number" data-line-number="631"></td>
        <td id="LC631" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzim;</td>
      </tr>
      <tr>
        <td id="L632" class="blob-num js-line-number" data-line-number="632"></td>
        <td id="LC632" class="blob-code js-file-line">    <span class="pl-k">int</span> iNeighborhood;</td>
      </tr>
      <tr>
        <td id="L633" class="blob-num js-line-number" data-line-number="633"></td>
        <td id="LC633" class="blob-code js-file-line">    <span class="pl-k">int</span> nNeighborhood;</td>
      </tr>
      <tr>
        <td id="L634" class="blob-num js-line-number" data-line-number="634"></td>
        <td id="LC634" class="blob-code js-file-line">    <span class="pl-k">int</span> count;</td>
      </tr>
      <tr>
        <td id="L635" class="blob-num js-line-number" data-line-number="635"></td>
        <td id="LC635" class="blob-code js-file-line">    <span class="pl-k">unsigned</span> <span class="pl-k">char</span> vradMissingValue;</td>
      </tr>
      <tr>
        <td id="L636" class="blob-num js-line-number" data-line-number="636"></td>
        <td id="LC636" class="blob-code js-file-line">    <span class="pl-k">unsigned</span> <span class="pl-k">char</span> dbzMissingValue;</td>
      </tr>
      <tr>
        <td id="L637" class="blob-num js-line-number" data-line-number="637"></td>
        <td id="LC637" class="blob-code js-file-line">    <span class="pl-k">unsigned</span> <span class="pl-k">char</span> texMissingValue;</td>
      </tr>
      <tr>
        <td id="L638" class="blob-num js-line-number" data-line-number="638"></td>
        <td id="LC638" class="blob-code js-file-line">    <span class="pl-k">double</span> vmoment1;</td>
      </tr>
      <tr>
        <td id="L639" class="blob-num js-line-number" data-line-number="639"></td>
        <td id="LC639" class="blob-code js-file-line">    <span class="pl-k">double</span> vmoment2;</td>
      </tr>
      <tr>
        <td id="L640" class="blob-num js-line-number" data-line-number="640"></td>
        <td id="LC640" class="blob-code js-file-line">    <span class="pl-k">double</span> dbz;</td>
      </tr>
      <tr>
        <td id="L641" class="blob-num js-line-number" data-line-number="641"></td>
        <td id="LC641" class="blob-code js-file-line">    <span class="pl-k">double</span> tex;</td>
      </tr>
      <tr>
        <td id="L642" class="blob-num js-line-number" data-line-number="642"></td>
        <td id="LC642" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L643" class="blob-num js-line-number" data-line-number="643"></td>
        <td id="LC643" class="blob-code js-file-line">    <span class="pl-k">int</span> iLocal;</td>
      </tr>
      <tr>
        <td id="L644" class="blob-num js-line-number" data-line-number="644"></td>
        <td id="LC644" class="blob-code js-file-line">    <span class="pl-k">float</span> texOffset;</td>
      </tr>
      <tr>
        <td id="L645" class="blob-num js-line-number" data-line-number="645"></td>
        <td id="LC645" class="blob-code js-file-line">    <span class="pl-k">float</span> texScale;</td>
      </tr>
      <tr>
        <td id="L646" class="blob-num js-line-number" data-line-number="646"></td>
        <td id="LC646" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzOffset;</td>
      </tr>
      <tr>
        <td id="L647" class="blob-num js-line-number" data-line-number="647"></td>
        <td id="LC647" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzScale;</td>
      </tr>
      <tr>
        <td id="L648" class="blob-num js-line-number" data-line-number="648"></td>
        <td id="LC648" class="blob-code js-file-line">    <span class="pl-k">float</span> vradOffset;</td>
      </tr>
      <tr>
        <td id="L649" class="blob-num js-line-number" data-line-number="649"></td>
        <td id="LC649" class="blob-code js-file-line">    <span class="pl-k">float</span> vradScale;</td>
      </tr>
      <tr>
        <td id="L650" class="blob-num js-line-number" data-line-number="650"></td>
        <td id="LC650" class="blob-code js-file-line">    <span class="pl-k">float</span> vRadDiff;</td>
      </tr>
      <tr>
        <td id="L651" class="blob-num js-line-number" data-line-number="651"></td>
        <td id="LC651" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L652" class="blob-num js-line-number" data-line-number="652"></td>
        <td id="LC652" class="blob-code js-file-line">    nRang = vradMeta-&gt;nRang;</td>
      </tr>
      <tr>
        <td id="L653" class="blob-num js-line-number" data-line-number="653"></td>
        <td id="LC653" class="blob-code js-file-line">    nAzim = vradMeta-&gt;nAzim;</td>
      </tr>
      <tr>
        <td id="L654" class="blob-num js-line-number" data-line-number="654"></td>
        <td id="LC654" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L655" class="blob-num js-line-number" data-line-number="655"></td>
        <td id="LC655" class="blob-code js-file-line">    dbzOffset = dbzMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L656" class="blob-num js-line-number" data-line-number="656"></td>
        <td id="LC656" class="blob-code js-file-line">    dbzScale = dbzMeta-&gt;valueScale;</td>
      </tr>
      <tr>
        <td id="L657" class="blob-num js-line-number" data-line-number="657"></td>
        <td id="LC657" class="blob-code js-file-line">    dbzMissingValue = dbzMeta-&gt;missing;</td>
      </tr>
      <tr>
        <td id="L658" class="blob-num js-line-number" data-line-number="658"></td>
        <td id="LC658" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L659" class="blob-num js-line-number" data-line-number="659"></td>
        <td id="LC659" class="blob-code js-file-line">    vradOffset = vradMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L660" class="blob-num js-line-number" data-line-number="660"></td>
        <td id="LC660" class="blob-code js-file-line">    vradScale = vradMeta-&gt;valueScale;</td>
      </tr>
      <tr>
        <td id="L661" class="blob-num js-line-number" data-line-number="661"></td>
        <td id="LC661" class="blob-code js-file-line">    vradMissingValue = vradMeta-&gt;missing;</td>
      </tr>
      <tr>
        <td id="L662" class="blob-num js-line-number" data-line-number="662"></td>
        <td id="LC662" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L663" class="blob-num js-line-number" data-line-number="663"></td>
        <td id="LC663" class="blob-code js-file-line">    texOffset = texMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L664" class="blob-num js-line-number" data-line-number="664"></td>
        <td id="LC664" class="blob-code js-file-line">    texScale = texMeta-&gt;valueScale;</td>
      </tr>
      <tr>
        <td id="L665" class="blob-num js-line-number" data-line-number="665"></td>
        <td id="LC665" class="blob-code js-file-line">    texMissingValue = texMeta-&gt;missing;</td>
      </tr>
      <tr>
        <td id="L666" class="blob-num js-line-number" data-line-number="666"></td>
        <td id="LC666" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L667" class="blob-num js-line-number" data-line-number="667"></td>
        <td id="LC667" class="blob-code js-file-line">    nNeighborhood = (<span class="pl-k">int</span>)(nRangNeighborhood * nAzimNeighborhood);</td>
      </tr>
      <tr>
        <td id="L668" class="blob-num js-line-number" data-line-number="668"></td>
        <td id="LC668" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L669" class="blob-num js-line-number" data-line-number="669"></td>
        <td id="LC669" class="blob-code js-file-line">    <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) {</td>
      </tr>
      <tr>
        <td id="L670" class="blob-num js-line-number" data-line-number="670"></td>
        <td id="LC670" class="blob-code js-file-line">        <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L671" class="blob-num js-line-number" data-line-number="671"></td>
        <td id="LC671" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L672" class="blob-num js-line-number" data-line-number="672"></td>
        <td id="LC672" class="blob-code js-file-line">            iGlobal = iRang + iAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L673" class="blob-num js-line-number" data-line-number="673"></td>
        <td id="LC673" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L674" class="blob-num js-line-number" data-line-number="674"></td>
        <td id="LC674" class="blob-code js-file-line">            <span class="pl-c">// count number of direct neighbors above threshold</span></td>
      </tr>
      <tr>
        <td id="L675" class="blob-num js-line-number" data-line-number="675"></td>
        <td id="LC675" class="blob-code js-file-line">            count = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L676" class="blob-num js-line-number" data-line-number="676"></td>
        <td id="LC676" class="blob-code js-file-line">            vmoment1 = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L677" class="blob-num js-line-number" data-line-number="677"></td>
        <td id="LC677" class="blob-code js-file-line">            vmoment2 = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L678" class="blob-num js-line-number" data-line-number="678"></td>
        <td id="LC678" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L679" class="blob-num js-line-number" data-line-number="679"></td>
        <td id="LC679" class="blob-code js-file-line">            dbz = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L680" class="blob-num js-line-number" data-line-number="680"></td>
        <td id="LC680" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L681" class="blob-num js-line-number" data-line-number="681"></td>
        <td id="LC681" class="blob-code js-file-line">            <span class="pl-k">for</span> (iNeighborhood = <span class="pl-c1">0</span>; iNeighborhood &lt; nNeighborhood; iNeighborhood++) {</td>
      </tr>
      <tr>
        <td id="L682" class="blob-num js-line-number" data-line-number="682"></td>
        <td id="LC682" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L683" class="blob-num js-line-number" data-line-number="683"></td>
        <td id="LC683" class="blob-code js-file-line">                iLocal = <span class="pl-c1">findNearbyGateIndex</span>(nAzim,nRang,iGlobal,nAzimNeighborhood,nRangNeighborhood,iNeighborhood);</td>
      </tr>
      <tr>
        <td id="L684" class="blob-num js-line-number" data-line-number="684"></td>
        <td id="LC684" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L685" class="blob-num js-line-number" data-line-number="685"></td>
        <td id="LC685" class="blob-code js-file-line">                #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L686" class="blob-num js-line-number" data-line-number="686"></td>
        <td id="LC686" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>iLocal = <span class="pl-c1">%d</span>; <span class="pl-pds">&quot;</span></span>,iLocal);</td>
      </tr>
      <tr>
        <td id="L687" class="blob-num js-line-number" data-line-number="687"></td>
        <td id="LC687" class="blob-code js-file-line">                #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L688" class="blob-num js-line-number" data-line-number="688"></td>
        <td id="LC688" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L689" class="blob-num js-line-number" data-line-number="689"></td>
        <td id="LC689" class="blob-code js-file-line">                <span class="pl-k">if</span> (iLocal &lt; <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L690" class="blob-num js-line-number" data-line-number="690"></td>
        <td id="LC690" class="blob-code js-file-line">                    <span class="pl-c">// iLocal less than zero are error codes</span></td>
      </tr>
      <tr>
        <td id="L691" class="blob-num js-line-number" data-line-number="691"></td>
        <td id="LC691" class="blob-code js-file-line">                    <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L692" class="blob-num js-line-number" data-line-number="692"></td>
        <td id="LC692" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L693" class="blob-num js-line-number" data-line-number="693"></td>
        <td id="LC693" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L694" class="blob-num js-line-number" data-line-number="694"></td>
        <td id="LC694" class="blob-code js-file-line">                <span class="pl-k">if</span> (vradImage[iLocal] == vradMissingValue || dbzImage[iLocal] == dbzMissingValue) {</td>
      </tr>
      <tr>
        <td id="L695" class="blob-num js-line-number" data-line-number="695"></td>
        <td id="LC695" class="blob-code js-file-line">                    <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L696" class="blob-num js-line-number" data-line-number="696"></td>
        <td id="LC696" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L697" class="blob-num js-line-number" data-line-number="697"></td>
        <td id="LC697" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L698" class="blob-num js-line-number" data-line-number="698"></td>
        <td id="LC698" class="blob-code js-file-line">                vRadDiff = vradOffset + vradScale * (<span class="pl-k">float</span>) (vradImage[iGlobal] - vradImage[iLocal]);</td>
      </tr>
      <tr>
        <td id="L699" class="blob-num js-line-number" data-line-number="699"></td>
        <td id="LC699" class="blob-code js-file-line">                vmoment1 += vRadDiff;</td>
      </tr>
      <tr>
        <td id="L700" class="blob-num js-line-number" data-line-number="700"></td>
        <td id="LC700" class="blob-code js-file-line">                vmoment2 += <span class="pl-c1">SQUARE</span>(vRadDiff);</td>
      </tr>
      <tr>
        <td id="L701" class="blob-num js-line-number" data-line-number="701"></td>
        <td id="LC701" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L702" class="blob-num js-line-number" data-line-number="702"></td>
        <td id="LC702" class="blob-code js-file-line">                dbz += dbzOffset + dbzScale * (<span class="pl-k">float</span>) dbzImage[iLocal];</td>
      </tr>
      <tr>
        <td id="L703" class="blob-num js-line-number" data-line-number="703"></td>
        <td id="LC703" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L704" class="blob-num js-line-number" data-line-number="704"></td>
        <td id="LC704" class="blob-code js-file-line">                count++;</td>
      </tr>
      <tr>
        <td id="L705" class="blob-num js-line-number" data-line-number="705"></td>
        <td id="LC705" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L706" class="blob-num js-line-number" data-line-number="706"></td>
        <td id="LC706" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L707" class="blob-num js-line-number" data-line-number="707"></td>
        <td id="LC707" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L708" class="blob-num js-line-number" data-line-number="708"></td>
        <td id="LC708" class="blob-code js-file-line">            vmoment1 /= count;</td>
      </tr>
      <tr>
        <td id="L709" class="blob-num js-line-number" data-line-number="709"></td>
        <td id="LC709" class="blob-code js-file-line">            vmoment2 /= count;</td>
      </tr>
      <tr>
        <td id="L710" class="blob-num js-line-number" data-line-number="710"></td>
        <td id="LC710" class="blob-code js-file-line">            dbz /= count;</td>
      </tr>
      <tr>
        <td id="L711" class="blob-num js-line-number" data-line-number="711"></td>
        <td id="LC711" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L712" class="blob-num js-line-number" data-line-number="712"></td>
        <td id="LC712" class="blob-code js-file-line">            <span class="pl-c">// when not enough neighbors, continue</span></td>
      </tr>
      <tr>
        <td id="L713" class="blob-num js-line-number" data-line-number="713"></td>
        <td id="LC713" class="blob-code js-file-line">            <span class="pl-k">if</span> (count &lt; nCountMin) {</td>
      </tr>
      <tr>
        <td id="L714" class="blob-num js-line-number" data-line-number="714"></td>
        <td id="LC714" class="blob-code js-file-line">                texImage[iGlobal] = texMissingValue;</td>
      </tr>
      <tr>
        <td id="L715" class="blob-num js-line-number" data-line-number="715"></td>
        <td id="LC715" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L716" class="blob-num js-line-number" data-line-number="716"></td>
        <td id="LC716" class="blob-code js-file-line">            <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L717" class="blob-num js-line-number" data-line-number="717"></td>
        <td id="LC717" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L718" class="blob-num js-line-number" data-line-number="718"></td>
        <td id="LC718" class="blob-code js-file-line">                tex = <span class="pl-c1">sqrt</span>(<span class="pl-c1">XABS</span>(vmoment2-<span class="pl-c1">SQUARE</span>(vmoment1)));</td>
      </tr>
      <tr>
        <td id="L719" class="blob-num js-line-number" data-line-number="719"></td>
        <td id="LC719" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L720" class="blob-num js-line-number" data-line-number="720"></td>
        <td id="LC720" class="blob-code js-file-line">                <span class="pl-k">double</span> tmpTex = <span class="pl-c1">ROUND</span>((tex - texOffset) / texScale);</td>
      </tr>
      <tr>
        <td id="L721" class="blob-num js-line-number" data-line-number="721"></td>
        <td id="LC721" class="blob-code js-file-line">                <span class="pl-k">if</span> (<span class="pl-c1">0</span> &lt;= tmpTex &amp;&amp; tmpTex &lt;= <span class="pl-c1">255</span>) {</td>
      </tr>
      <tr>
        <td id="L722" class="blob-num js-line-number" data-line-number="722"></td>
        <td id="LC722" class="blob-code js-file-line">                    texImage[iGlobal] = (<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) tmpTex;</td>
      </tr>
      <tr>
        <td id="L723" class="blob-num js-line-number" data-line-number="723"></td>
        <td id="LC723" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L724" class="blob-num js-line-number" data-line-number="724"></td>
        <td id="LC724" class="blob-code js-file-line">                <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L725" class="blob-num js-line-number" data-line-number="725"></td>
        <td id="LC725" class="blob-code js-file-line">                    <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Error casting texture value of <span class="pl-c1">%f</span> to unsigned char type at texImage[<span class="pl-c1">%d</span>]. Aborting.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,tmpTex,iGlobal);</td>
      </tr>
      <tr>
        <td id="L726" class="blob-num js-line-number" data-line-number="726"></td>
        <td id="LC726" class="blob-code js-file-line">                    <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L727" class="blob-num js-line-number" data-line-number="727"></td>
        <td id="LC727" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L728" class="blob-num js-line-number" data-line-number="728"></td>
        <td id="LC728" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L729" class="blob-num js-line-number" data-line-number="729"></td>
        <td id="LC729" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L730" class="blob-num js-line-number" data-line-number="730"></td>
        <td id="LC730" class="blob-code js-file-line">                #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L731" class="blob-num js-line-number" data-line-number="731"></td>
        <td id="LC731" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr,</td>
      </tr>
      <tr>
        <td id="L732" class="blob-num js-line-number" data-line-number="732"></td>
        <td id="LC732" class="blob-code js-file-line">                        <span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n</span>(C) count = <span class="pl-c1">%d</span>; nCountMin = <span class="pl-c1">%d</span>; vmoment1 = <span class="pl-c1">%f</span>; vmoment2 = <span class="pl-c1">%f</span>; tex = <span class="pl-c1">%f</span>; texBody[<span class="pl-c1">%d</span>] = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L733" class="blob-num js-line-number" data-line-number="733"></td>
        <td id="LC733" class="blob-code js-file-line">                        count, nCountMin, vmoment1, vmoment2, tex,</td>
      </tr>
      <tr>
        <td id="L734" class="blob-num js-line-number" data-line-number="734"></td>
        <td id="LC734" class="blob-code js-file-line">                        iGlobal, texImage[iGlobal]);</td>
      </tr>
      <tr>
        <td id="L735" class="blob-num js-line-number" data-line-number="735"></td>
        <td id="LC735" class="blob-code js-file-line">                #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L736" class="blob-num js-line-number" data-line-number="736"></td>
        <td id="LC736" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L737" class="blob-num js-line-number" data-line-number="737"></td>
        <td id="LC737" class="blob-code js-file-line">            } <span class="pl-c">//else</span></td>
      </tr>
      <tr>
        <td id="L738" class="blob-num js-line-number" data-line-number="738"></td>
        <td id="LC738" class="blob-code js-file-line">        } <span class="pl-c">//for</span></td>
      </tr>
      <tr>
        <td id="L739" class="blob-num js-line-number" data-line-number="739"></td>
        <td id="LC739" class="blob-code js-file-line">    } <span class="pl-c">//for</span></td>
      </tr>
      <tr>
        <td id="L740" class="blob-num js-line-number" data-line-number="740"></td>
        <td id="LC740" class="blob-code js-file-line">} <span class="pl-c">// calcTexture</span></td>
      </tr>
      <tr>
        <td id="L741" class="blob-num js-line-number" data-line-number="741"></td>
        <td id="LC741" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L742" class="blob-num js-line-number" data-line-number="742"></td>
        <td id="LC742" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L743" class="blob-num js-line-number" data-line-number="743"></td>
        <td id="LC743" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L744" class="blob-num js-line-number" data-line-number="744"></td>
        <td id="LC744" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L745" class="blob-num js-line-number" data-line-number="745"></td>
        <td id="LC745" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L746" class="blob-num js-line-number" data-line-number="746"></td>
        <td id="LC746" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">classifyGatesSimple</span>(<span class="pl-k">void</span>) {</td>
      </tr>
      <tr>
        <td id="L747" class="blob-num js-line-number" data-line-number="747"></td>
        <td id="LC747" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L748" class="blob-num js-line-number" data-line-number="748"></td>
        <td id="LC748" class="blob-code js-file-line">    <span class="pl-k">int</span> iPoint;</td>
      </tr>
      <tr>
        <td id="L749" class="blob-num js-line-number" data-line-number="749"></td>
        <td id="LC749" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L750" class="blob-num js-line-number" data-line-number="750"></td>
        <td id="LC750" class="blob-code js-file-line">    <span class="pl-k">for</span> (iPoint = <span class="pl-c1">0</span>; iPoint &lt; nRowsPoints; iPoint++) {</td>
      </tr>
      <tr>
        <td id="L751" class="blob-num js-line-number" data-line-number="751"></td>
        <td id="LC751" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L752" class="blob-num js-line-number" data-line-number="752"></td>
        <td id="LC752" class="blob-code js-file-line">        <span class="pl-k">const</span> <span class="pl-k">float</span> azimValue = points[iPoint * nColsPoints + azimAngleCol];</td>
      </tr>
      <tr>
        <td id="L753" class="blob-num js-line-number" data-line-number="753"></td>
        <td id="LC753" class="blob-code js-file-line">        <span class="pl-k">const</span> <span class="pl-k">float</span> dbzValue = points[iPoint * nColsPoints + dbzValueCol];        </td>
      </tr>
      <tr>
        <td id="L754" class="blob-num js-line-number" data-line-number="754"></td>
        <td id="LC754" class="blob-code js-file-line">        <span class="pl-k">const</span> <span class="pl-k">float</span> vradValue = points[iPoint * nColsPoints + vradValueCol];</td>
      </tr>
      <tr>
        <td id="L755" class="blob-num js-line-number" data-line-number="755"></td>
        <td id="LC755" class="blob-code js-file-line">        <span class="pl-k">const</span> <span class="pl-k">int</span> cellValue = (<span class="pl-k">int</span>) points[iPoint * nColsPoints + cellValueCol];</td>
      </tr>
      <tr>
        <td id="L756" class="blob-num js-line-number" data-line-number="756"></td>
        <td id="LC756" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L757" class="blob-num js-line-number" data-line-number="757"></td>
        <td id="LC757" class="blob-code js-file-line">        <span class="pl-k">unsigned</span> <span class="pl-k">int</span> gateCode = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L758" class="blob-num js-line-number" data-line-number="758"></td>
        <td id="LC758" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L759" class="blob-num js-line-number" data-line-number="759"></td>
        <td id="LC759" class="blob-code js-file-line">        <span class="pl-k">if</span> (<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L760" class="blob-num js-line-number" data-line-number="760"></td>
        <td id="LC760" class="blob-code js-file-line">            <span class="pl-c">// this gate is true in the static clutter map (which we don&#39;t have yet TODO)</span></td>
      </tr>
      <tr>
        <td id="L761" class="blob-num js-line-number" data-line-number="761"></td>
        <td id="LC761" class="blob-code js-file-line">            gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionStaticClutter;</td>
      </tr>
      <tr>
        <td id="L762" class="blob-num js-line-number" data-line-number="762"></td>
        <td id="LC762" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L763" class="blob-num js-line-number" data-line-number="763"></td>
        <td id="LC763" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L764" class="blob-num js-line-number" data-line-number="764"></td>
        <td id="LC764" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellValue == <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L765" class="blob-num js-line-number" data-line-number="765"></td>
        <td id="LC765" class="blob-code js-file-line">            <span class="pl-c">// this gate is part of the cluttermap (without fringe)</span></td>
      </tr>
      <tr>
        <td id="L766" class="blob-num js-line-number" data-line-number="766"></td>
        <td id="LC766" class="blob-code js-file-line">            gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionDynamicClutter;</td>
      </tr>
      <tr>
        <td id="L767" class="blob-num js-line-number" data-line-number="767"></td>
        <td id="LC767" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L768" class="blob-num js-line-number" data-line-number="768"></td>
        <td id="LC768" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L769" class="blob-num js-line-number" data-line-number="769"></td>
        <td id="LC769" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellValue == <span class="pl-c1">2</span>) {</td>
      </tr>
      <tr>
        <td id="L770" class="blob-num js-line-number" data-line-number="770"></td>
        <td id="LC770" class="blob-code js-file-line">            <span class="pl-c">// this gate is part of the fringe of the cluttermap</span></td>
      </tr>
      <tr>
        <td id="L771" class="blob-num js-line-number" data-line-number="771"></td>
        <td id="LC771" class="blob-code js-file-line">            gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionDynamicClutterFringe;</td>
      </tr>
      <tr>
        <td id="L772" class="blob-num js-line-number" data-line-number="772"></td>
        <td id="LC772" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L773" class="blob-num js-line-number" data-line-number="773"></td>
        <td id="LC773" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L774" class="blob-num js-line-number" data-line-number="774"></td>
        <td id="LC774" class="blob-code js-file-line">        <span class="pl-k">if</span> (<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L775" class="blob-num js-line-number" data-line-number="775"></td>
        <td id="LC775" class="blob-code js-file-line">            <span class="pl-c">// this gate has reflectivity data but no corresponding radial velocity data</span></td>
      </tr>
      <tr>
        <td id="L776" class="blob-num js-line-number" data-line-number="776"></td>
        <td id="LC776" class="blob-code js-file-line">            <span class="pl-c">// TODO no condition for this yet</span></td>
      </tr>
      <tr>
        <td id="L777" class="blob-num js-line-number" data-line-number="777"></td>
        <td id="LC777" class="blob-code js-file-line">            gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionVradMissing;</td>
      </tr>
      <tr>
        <td id="L778" class="blob-num js-line-number" data-line-number="778"></td>
        <td id="LC778" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L779" class="blob-num js-line-number" data-line-number="779"></td>
        <td id="LC779" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L780" class="blob-num js-line-number" data-line-number="780"></td>
        <td id="LC780" class="blob-code js-file-line">        <span class="pl-k">if</span> (dbzValue &gt; dbzMax) {</td>
      </tr>
      <tr>
        <td id="L781" class="blob-num js-line-number" data-line-number="781"></td>
        <td id="LC781" class="blob-code js-file-line">            <span class="pl-c">// this gate&#39;s dbz value is too high to be due to birds, it must be </span></td>
      </tr>
      <tr>
        <td id="L782" class="blob-num js-line-number" data-line-number="782"></td>
        <td id="LC782" class="blob-code js-file-line">            <span class="pl-c">// caused by something else</span></td>
      </tr>
      <tr>
        <td id="L783" class="blob-num js-line-number" data-line-number="783"></td>
        <td id="LC783" class="blob-code js-file-line">            gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionDbzTooHighForBirds;</td>
      </tr>
      <tr>
        <td id="L784" class="blob-num js-line-number" data-line-number="784"></td>
        <td id="LC784" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L785" class="blob-num js-line-number" data-line-number="785"></td>
        <td id="LC785" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L786" class="blob-num js-line-number" data-line-number="786"></td>
        <td id="LC786" class="blob-code js-file-line">        <span class="pl-k">if</span> (vradValue &lt; vradMin) {</td>
      </tr>
      <tr>
        <td id="L787" class="blob-num js-line-number" data-line-number="787"></td>
        <td id="LC787" class="blob-code js-file-line">            <span class="pl-c">// this gate&#39;s radial velocity is too low to be due to actual scatterers; likely just noise</span></td>
      </tr>
      <tr>
        <td id="L788" class="blob-num js-line-number" data-line-number="788"></td>
        <td id="LC788" class="blob-code js-file-line">            gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionVradTooLow;</td>
      </tr>
      <tr>
        <td id="L789" class="blob-num js-line-number" data-line-number="789"></td>
        <td id="LC789" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L790" class="blob-num js-line-number" data-line-number="790"></td>
        <td id="LC790" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L791" class="blob-num js-line-number" data-line-number="791"></td>
        <td id="LC791" class="blob-code js-file-line">        <span class="pl-k">if</span> (<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L792" class="blob-num js-line-number" data-line-number="792"></td>
        <td id="LC792" class="blob-code js-file-line">            <span class="pl-c">// this flag is set later on by updateFlagFieldsInPointsArray()</span></td>
      </tr>
      <tr>
        <td id="L793" class="blob-num js-line-number" data-line-number="793"></td>
        <td id="LC793" class="blob-code js-file-line">            <span class="pl-c">// flagPositionVDifMax</span></td>
      </tr>
      <tr>
        <td id="L794" class="blob-num js-line-number" data-line-number="794"></td>
        <td id="LC794" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L795" class="blob-num js-line-number" data-line-number="795"></td>
        <td id="LC795" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L796" class="blob-num js-line-number" data-line-number="796"></td>
        <td id="LC796" class="blob-code js-file-line">        <span class="pl-k">if</span> (azimValue &lt; azimMin) {</td>
      </tr>
      <tr>
        <td id="L797" class="blob-num js-line-number" data-line-number="797"></td>
        <td id="LC797" class="blob-code js-file-line">            <span class="pl-c">// the user can specify to exclude gates based on their azimuth;</span></td>
      </tr>
      <tr>
        <td id="L798" class="blob-num js-line-number" data-line-number="798"></td>
        <td id="LC798" class="blob-code js-file-line">            <span class="pl-c">// this clause is for gates that have too low azimuth</span></td>
      </tr>
      <tr>
        <td id="L799" class="blob-num js-line-number" data-line-number="799"></td>
        <td id="LC799" class="blob-code js-file-line">            gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionAzimTooLow;</td>
      </tr>
      <tr>
        <td id="L800" class="blob-num js-line-number" data-line-number="800"></td>
        <td id="LC800" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L801" class="blob-num js-line-number" data-line-number="801"></td>
        <td id="LC801" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L802" class="blob-num js-line-number" data-line-number="802"></td>
        <td id="LC802" class="blob-code js-file-line">        <span class="pl-k">if</span> (azimValue &gt; azimMax) {</td>
      </tr>
      <tr>
        <td id="L803" class="blob-num js-line-number" data-line-number="803"></td>
        <td id="LC803" class="blob-code js-file-line">            <span class="pl-c">// the user can specify to exclude gates based on their azimuth;</span></td>
      </tr>
      <tr>
        <td id="L804" class="blob-num js-line-number" data-line-number="804"></td>
        <td id="LC804" class="blob-code js-file-line">            <span class="pl-c">// this clause is for gates that have too high azimuth</span></td>
      </tr>
      <tr>
        <td id="L805" class="blob-num js-line-number" data-line-number="805"></td>
        <td id="LC805" class="blob-code js-file-line">            gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionAzimTooHigh;</td>
      </tr>
      <tr>
        <td id="L806" class="blob-num js-line-number" data-line-number="806"></td>
        <td id="LC806" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L807" class="blob-num js-line-number" data-line-number="807"></td>
        <td id="LC807" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L808" class="blob-num js-line-number" data-line-number="808"></td>
        <td id="LC808" class="blob-code js-file-line">        points[iPoint * nColsPoints + gateCodeCol] = (<span class="pl-k">float</span>) gateCode;</td>
      </tr>
      <tr>
        <td id="L809" class="blob-num js-line-number" data-line-number="809"></td>
        <td id="LC809" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L810" class="blob-num js-line-number" data-line-number="810"></td>
        <td id="LC810" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L811" class="blob-num js-line-number" data-line-number="811"></td>
        <td id="LC811" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L812" class="blob-num js-line-number" data-line-number="812"></td>
        <td id="LC812" class="blob-code js-file-line">    <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L813" class="blob-num js-line-number" data-line-number="813"></td>
        <td id="LC813" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L814" class="blob-num js-line-number" data-line-number="814"></td>
        <td id="LC814" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L815" class="blob-num js-line-number" data-line-number="815"></td>
        <td id="LC815" class="blob-code js-file-line">};</td>
      </tr>
      <tr>
        <td id="L816" class="blob-num js-line-number" data-line-number="816"></td>
        <td id="LC816" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L817" class="blob-num js-line-number" data-line-number="817"></td>
        <td id="LC817" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L818" class="blob-num js-line-number" data-line-number="818"></td>
        <td id="LC818" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L819" class="blob-num js-line-number" data-line-number="819"></td>
        <td id="LC819" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L820" class="blob-num js-line-number" data-line-number="820"></td>
        <td id="LC820" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L821" class="blob-num js-line-number" data-line-number="821"></td>
        <td id="LC821" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">constructorInt</span>(SCANMETA* meta, <span class="pl-k">int</span>* image, PolarScan_t* scan, <span class="pl-k">const</span> <span class="pl-k">int</span> nGlobal, <span class="pl-k">const</span> <span class="pl-k">int</span> initValue) {</td>
      </tr>
      <tr>
        <td id="L822" class="blob-num js-line-number" data-line-number="822"></td>
        <td id="LC822" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L823" class="blob-num js-line-number" data-line-number="823"></td>
        <td id="LC823" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L824" class="blob-num js-line-number" data-line-number="824"></td>
        <td id="LC824" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L825" class="blob-num js-line-number" data-line-number="825"></td>
        <td id="LC825" class="blob-code js-file-line">    <span class="pl-k">if</span> (image == <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L826" class="blob-num js-line-number" data-line-number="826"></td>
        <td id="LC826" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Error allocating image<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L827" class="blob-num js-line-number" data-line-number="827"></td>
        <td id="LC827" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">2</span>;</td>
      </tr>
      <tr>
        <td id="L828" class="blob-num js-line-number" data-line-number="828"></td>
        <td id="LC828" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L829" class="blob-num js-line-number" data-line-number="829"></td>
        <td id="LC829" class="blob-code js-file-line">    <span class="pl-k">for</span> (iGlobal = <span class="pl-c1">0</span>; iGlobal &lt; nGlobal; iGlobal++) {</td>
      </tr>
      <tr>
        <td id="L830" class="blob-num js-line-number" data-line-number="830"></td>
        <td id="LC830" class="blob-code js-file-line">        image[iGlobal] = initValue;</td>
      </tr>
      <tr>
        <td id="L831" class="blob-num js-line-number" data-line-number="831"></td>
        <td id="LC831" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L832" class="blob-num js-line-number" data-line-number="832"></td>
        <td id="LC832" class="blob-code js-file-line">    meta-&gt;heig = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScan_getHeight</span>(scan);</td>
      </tr>
      <tr>
        <td id="L833" class="blob-num js-line-number" data-line-number="833"></td>
        <td id="LC833" class="blob-code js-file-line">    meta-&gt;elev = (<span class="pl-k">float</span>) (<span class="pl-c1">360</span> * <span class="pl-c1">PolarScan_getElangle</span>(scan) / <span class="pl-c1">2</span> / PI );</td>
      </tr>
      <tr>
        <td id="L834" class="blob-num js-line-number" data-line-number="834"></td>
        <td id="LC834" class="blob-code js-file-line">    meta-&gt;nRang = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNbins</span>(scan);</td>
      </tr>
      <tr>
        <td id="L835" class="blob-num js-line-number" data-line-number="835"></td>
        <td id="LC835" class="blob-code js-file-line">    meta-&gt;nAzim = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNrays</span>(scan);</td>
      </tr>
      <tr>
        <td id="L836" class="blob-num js-line-number" data-line-number="836"></td>
        <td id="LC836" class="blob-code js-file-line">    meta-&gt;rangeScale = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScan_getRscale</span>(scan);</td>
      </tr>
      <tr>
        <td id="L837" class="blob-num js-line-number" data-line-number="837"></td>
        <td id="LC837" class="blob-code js-file-line">    meta-&gt;azimScale = <span class="pl-c1">360</span>.<span class="pl-c1">0f</span>/meta-&gt;nAzim;   <span class="pl-c">// FIXME not sure if this always works</span></td>
      </tr>
      <tr>
        <td id="L838" class="blob-num js-line-number" data-line-number="838"></td>
        <td id="LC838" class="blob-code js-file-line">    meta-&gt;valueOffset = <span class="pl-c1">0</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L839" class="blob-num js-line-number" data-line-number="839"></td>
        <td id="LC839" class="blob-code js-file-line">    meta-&gt;valueScale = <span class="pl-c1">1</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L840" class="blob-num js-line-number" data-line-number="840"></td>
        <td id="LC840" class="blob-code js-file-line">    meta-&gt;missing = (<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) <span class="pl-c1">255</span>;  <span class="pl-c">// FIXME this does not work as intended for type int</span></td>
      </tr>
      <tr>
        <td id="L841" class="blob-num js-line-number" data-line-number="841"></td>
        <td id="LC841" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L842" class="blob-num js-line-number" data-line-number="842"></td>
        <td id="LC842" class="blob-code js-file-line">    <span class="pl-k">return</span> <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L843" class="blob-num js-line-number" data-line-number="843"></td>
        <td id="LC843" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L844" class="blob-num js-line-number" data-line-number="844"></td>
        <td id="LC844" class="blob-code js-file-line">}</td>
      </tr>
      <tr>
        <td id="L845" class="blob-num js-line-number" data-line-number="845"></td>
        <td id="LC845" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L846" class="blob-num js-line-number" data-line-number="846"></td>
        <td id="LC846" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L847" class="blob-num js-line-number" data-line-number="847"></td>
        <td id="LC847" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L848" class="blob-num js-line-number" data-line-number="848"></td>
        <td id="LC848" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L849" class="blob-num js-line-number" data-line-number="849"></td>
        <td id="LC849" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">constructorUChar</span>(SCANMETA* meta, <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* image, PolarScan_t* scan, <span class="pl-k">const</span> <span class="pl-k">int</span> nGlobal, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> initValue) {</td>
      </tr>
      <tr>
        <td id="L850" class="blob-num js-line-number" data-line-number="850"></td>
        <td id="LC850" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L851" class="blob-num js-line-number" data-line-number="851"></td>
        <td id="LC851" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L852" class="blob-num js-line-number" data-line-number="852"></td>
        <td id="LC852" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L853" class="blob-num js-line-number" data-line-number="853"></td>
        <td id="LC853" class="blob-code js-file-line">    <span class="pl-k">if</span> (image == <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L854" class="blob-num js-line-number" data-line-number="854"></td>
        <td id="LC854" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Error allocating image<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L855" class="blob-num js-line-number" data-line-number="855"></td>
        <td id="LC855" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">2</span>;</td>
      </tr>
      <tr>
        <td id="L856" class="blob-num js-line-number" data-line-number="856"></td>
        <td id="LC856" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L857" class="blob-num js-line-number" data-line-number="857"></td>
        <td id="LC857" class="blob-code js-file-line">    <span class="pl-k">for</span> (iGlobal = <span class="pl-c1">0</span>; iGlobal &lt; nGlobal; iGlobal++) {</td>
      </tr>
      <tr>
        <td id="L858" class="blob-num js-line-number" data-line-number="858"></td>
        <td id="LC858" class="blob-code js-file-line">        image[iGlobal] = initValue;</td>
      </tr>
      <tr>
        <td id="L859" class="blob-num js-line-number" data-line-number="859"></td>
        <td id="LC859" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L860" class="blob-num js-line-number" data-line-number="860"></td>
        <td id="LC860" class="blob-code js-file-line">    meta-&gt;heig = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScan_getHeight</span>(scan);</td>
      </tr>
      <tr>
        <td id="L861" class="blob-num js-line-number" data-line-number="861"></td>
        <td id="LC861" class="blob-code js-file-line">    meta-&gt;elev = (<span class="pl-k">float</span>) (<span class="pl-c1">360</span> * <span class="pl-c1">PolarScan_getElangle</span>(scan) / <span class="pl-c1">2</span> / PI );</td>
      </tr>
      <tr>
        <td id="L862" class="blob-num js-line-number" data-line-number="862"></td>
        <td id="LC862" class="blob-code js-file-line">    meta-&gt;nRang = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNbins</span>(scan);</td>
      </tr>
      <tr>
        <td id="L863" class="blob-num js-line-number" data-line-number="863"></td>
        <td id="LC863" class="blob-code js-file-line">    meta-&gt;nAzim = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNrays</span>(scan);</td>
      </tr>
      <tr>
        <td id="L864" class="blob-num js-line-number" data-line-number="864"></td>
        <td id="LC864" class="blob-code js-file-line">    meta-&gt;rangeScale = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScan_getRscale</span>(scan);</td>
      </tr>
      <tr>
        <td id="L865" class="blob-num js-line-number" data-line-number="865"></td>
        <td id="LC865" class="blob-code js-file-line">    meta-&gt;azimScale = <span class="pl-c1">360</span>.<span class="pl-c1">0f</span>/meta-&gt;nAzim;   <span class="pl-c">// FIXME not sure if this always works</span></td>
      </tr>
      <tr>
        <td id="L866" class="blob-num js-line-number" data-line-number="866"></td>
        <td id="LC866" class="blob-code js-file-line">    meta-&gt;valueOffset = <span class="pl-c1">0</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L867" class="blob-num js-line-number" data-line-number="867"></td>
        <td id="LC867" class="blob-code js-file-line">    meta-&gt;valueScale = <span class="pl-c1">1</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L868" class="blob-num js-line-number" data-line-number="868"></td>
        <td id="LC868" class="blob-code js-file-line">    meta-&gt;missing = (<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) <span class="pl-c1">255</span>;</td>
      </tr>
      <tr>
        <td id="L869" class="blob-num js-line-number" data-line-number="869"></td>
        <td id="LC869" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L870" class="blob-num js-line-number" data-line-number="870"></td>
        <td id="LC870" class="blob-code js-file-line">    <span class="pl-k">return</span> <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L871" class="blob-num js-line-number" data-line-number="871"></td>
        <td id="LC871" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L872" class="blob-num js-line-number" data-line-number="872"></td>
        <td id="LC872" class="blob-code js-file-line">}</td>
      </tr>
      <tr>
        <td id="L873" class="blob-num js-line-number" data-line-number="873"></td>
        <td id="LC873" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L874" class="blob-num js-line-number" data-line-number="874"></td>
        <td id="LC874" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L875" class="blob-num js-line-number" data-line-number="875"></td>
        <td id="LC875" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L876" class="blob-num js-line-number" data-line-number="876"></td>
        <td id="LC876" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">constructPointsArray</span>(PolarVolume_t* volume) {</td>
      </tr>
      <tr>
        <td id="L877" class="blob-num js-line-number" data-line-number="877"></td>
        <td id="LC877" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L878" class="blob-num js-line-number" data-line-number="878"></td>
        <td id="LC878" class="blob-code js-file-line">        <span class="pl-c">// iterate over the scans in &#39;volume&#39;</span></td>
      </tr>
      <tr>
        <td id="L879" class="blob-num js-line-number" data-line-number="879"></td>
        <td id="LC879" class="blob-code js-file-line">        <span class="pl-k">int</span> iScan;</td>
      </tr>
      <tr>
        <td id="L880" class="blob-num js-line-number" data-line-number="880"></td>
        <td id="LC880" class="blob-code js-file-line">        <span class="pl-k">int</span> nScans;</td>
      </tr>
      <tr>
        <td id="L881" class="blob-num js-line-number" data-line-number="881"></td>
        <td id="LC881" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L882" class="blob-num js-line-number" data-line-number="882"></td>
        <td id="LC882" class="blob-code js-file-line">        <span class="pl-c">// determine how many scan elevations the volume object contains</span></td>
      </tr>
      <tr>
        <td id="L883" class="blob-num js-line-number" data-line-number="883"></td>
        <td id="LC883" class="blob-code js-file-line">        nScans = <span class="pl-c1">PolarVolume_getNumberOfScans</span>(volume);</td>
      </tr>
      <tr>
        <td id="L884" class="blob-num js-line-number" data-line-number="884"></td>
        <td id="LC884" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L885" class="blob-num js-line-number" data-line-number="885"></td>
        <td id="LC885" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L886" class="blob-num js-line-number" data-line-number="886"></td>
        <td id="LC886" class="blob-code js-file-line">        <span class="pl-k">for</span> (iScan = <span class="pl-c1">0</span>; iScan &lt; nScans; iScan++) {</td>
      </tr>
      <tr>
        <td id="L887" class="blob-num js-line-number" data-line-number="887"></td>
        <td id="LC887" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L888" class="blob-num js-line-number" data-line-number="888"></td>
        <td id="LC888" class="blob-code js-file-line">            <span class="pl-c">// initialize the scan object</span></td>
      </tr>
      <tr>
        <td id="L889" class="blob-num js-line-number" data-line-number="889"></td>
        <td id="LC889" class="blob-code js-file-line">            PolarScan_t* scan = <span class="pl-c1">NULL</span>;</td>
      </tr>
      <tr>
        <td id="L890" class="blob-num js-line-number" data-line-number="890"></td>
        <td id="LC890" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L891" class="blob-num js-line-number" data-line-number="891"></td>
        <td id="LC891" class="blob-code js-file-line">            <span class="pl-c">// extract the scan object from the volume object</span></td>
      </tr>
      <tr>
        <td id="L892" class="blob-num js-line-number" data-line-number="892"></td>
        <td id="LC892" class="blob-code js-file-line">            scan = <span class="pl-c1">PolarVolume_getScan</span>(volume, iScan);</td>
      </tr>
      <tr>
        <td id="L893" class="blob-num js-line-number" data-line-number="893"></td>
        <td id="LC893" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L894" class="blob-num js-line-number" data-line-number="894"></td>
        <td id="LC894" class="blob-code js-file-line">            <span class="pl-c">// determine the number of array elements in the polar scan</span></td>
      </tr>
      <tr>
        <td id="L895" class="blob-num js-line-number" data-line-number="895"></td>
        <td id="LC895" class="blob-code js-file-line">            <span class="pl-k">int</span> nGlobal = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNbins</span>(scan) * <span class="pl-c1">PolarScan_getNrays</span>(scan);</td>
      </tr>
      <tr>
        <td id="L896" class="blob-num js-line-number" data-line-number="896"></td>
        <td id="LC896" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L897" class="blob-num js-line-number" data-line-number="897"></td>
        <td id="LC897" class="blob-code js-file-line">            <span class="pl-c">// pre-allocate the dbz variables</span></td>
      </tr>
      <tr>
        <td id="L898" class="blob-num js-line-number" data-line-number="898"></td>
        <td id="LC898" class="blob-code js-file-line">            <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* dbzImage = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) * nGlobal);</td>
      </tr>
      <tr>
        <td id="L899" class="blob-num js-line-number" data-line-number="899"></td>
        <td id="LC899" class="blob-code js-file-line">            SCANMETA dbzMeta;</td>
      </tr>
      <tr>
        <td id="L900" class="blob-num js-line-number" data-line-number="900"></td>
        <td id="LC900" class="blob-code js-file-line">            <span class="pl-c1">constructorUChar</span>(&amp;dbzMeta, &amp;dbzImage[<span class="pl-c1">0</span>], scan, nGlobal, <span class="pl-c1">0</span>);</td>
      </tr>
      <tr>
        <td id="L901" class="blob-num js-line-number" data-line-number="901"></td>
        <td id="LC901" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L902" class="blob-num js-line-number" data-line-number="902"></td>
        <td id="LC902" class="blob-code js-file-line">            <span class="pl-c">// pre-allocate the vrad variables</span></td>
      </tr>
      <tr>
        <td id="L903" class="blob-num js-line-number" data-line-number="903"></td>
        <td id="LC903" class="blob-code js-file-line">            <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* vradImage = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) * nGlobal);</td>
      </tr>
      <tr>
        <td id="L904" class="blob-num js-line-number" data-line-number="904"></td>
        <td id="LC904" class="blob-code js-file-line">            SCANMETA vradMeta;</td>
      </tr>
      <tr>
        <td id="L905" class="blob-num js-line-number" data-line-number="905"></td>
        <td id="LC905" class="blob-code js-file-line">            <span class="pl-c1">constructorUChar</span>(&amp;vradMeta, &amp;vradImage[<span class="pl-c1">0</span>], scan, nGlobal, <span class="pl-c1">0</span>);</td>
      </tr>
      <tr>
        <td id="L906" class="blob-num js-line-number" data-line-number="906"></td>
        <td id="LC906" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L907" class="blob-num js-line-number" data-line-number="907"></td>
        <td id="LC907" class="blob-code js-file-line">            <span class="pl-c">// pre-allocate the tex variables</span></td>
      </tr>
      <tr>
        <td id="L908" class="blob-num js-line-number" data-line-number="908"></td>
        <td id="LC908" class="blob-code js-file-line">            <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* texImage = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) * nGlobal);</td>
      </tr>
      <tr>
        <td id="L909" class="blob-num js-line-number" data-line-number="909"></td>
        <td id="LC909" class="blob-code js-file-line">            SCANMETA texMeta;</td>
      </tr>
      <tr>
        <td id="L910" class="blob-num js-line-number" data-line-number="910"></td>
        <td id="LC910" class="blob-code js-file-line">            <span class="pl-c1">constructorUChar</span>(&amp;texMeta, &amp;texImage[<span class="pl-c1">0</span>], scan, nGlobal, <span class="pl-c1">0</span>);</td>
      </tr>
      <tr>
        <td id="L911" class="blob-num js-line-number" data-line-number="911"></td>
        <td id="LC911" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L912" class="blob-num js-line-number" data-line-number="912"></td>
        <td id="LC912" class="blob-code js-file-line">            <span class="pl-c">// pre-allocate the clutter variables</span></td>
      </tr>
      <tr>
        <td id="L913" class="blob-num js-line-number" data-line-number="913"></td>
        <td id="LC913" class="blob-code js-file-line">            <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* clutterImage = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) * nGlobal);</td>
      </tr>
      <tr>
        <td id="L914" class="blob-num js-line-number" data-line-number="914"></td>
        <td id="LC914" class="blob-code js-file-line">            SCANMETA clutterMeta;</td>
      </tr>
      <tr>
        <td id="L915" class="blob-num js-line-number" data-line-number="915"></td>
        <td id="LC915" class="blob-code js-file-line">            <span class="pl-c1">constructorUChar</span>(&amp;clutterMeta, &amp;clutterImage[<span class="pl-c1">0</span>], scan, nGlobal, <span class="pl-c1">0</span>);</td>
      </tr>
      <tr>
        <td id="L916" class="blob-num js-line-number" data-line-number="916"></td>
        <td id="LC916" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L917" class="blob-num js-line-number" data-line-number="917"></td>
        <td id="LC917" class="blob-code js-file-line">            <span class="pl-c">// pre-allocate the cell variables</span></td>
      </tr>
      <tr>
        <td id="L918" class="blob-num js-line-number" data-line-number="918"></td>
        <td id="LC918" class="blob-code js-file-line">            <span class="pl-k">int</span>* cellImage = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nGlobal);</td>
      </tr>
      <tr>
        <td id="L919" class="blob-num js-line-number" data-line-number="919"></td>
        <td id="LC919" class="blob-code js-file-line">            SCANMETA cellMeta;</td>
      </tr>
      <tr>
        <td id="L920" class="blob-num js-line-number" data-line-number="920"></td>
        <td id="LC920" class="blob-code js-file-line">            <span class="pl-c1">constructorInt</span>(&amp;cellMeta, &amp;cellImage[<span class="pl-c1">0</span>], scan, nGlobal, -<span class="pl-c1">1</span>);</td>
      </tr>
      <tr>
        <td id="L921" class="blob-num js-line-number" data-line-number="921"></td>
        <td id="LC921" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L922" class="blob-num js-line-number" data-line-number="922"></td>
        <td id="LC922" class="blob-code js-file-line">            <span class="pl-c">// populate the dbzMeta and dbzImage variables with data from </span></td>
      </tr>
      <tr>
        <td id="L923" class="blob-num js-line-number" data-line-number="923"></td>
        <td id="LC923" class="blob-code js-file-line">            <span class="pl-c">// the Rave scan object:</span></td>
      </tr>
      <tr>
        <td id="L924" class="blob-num js-line-number" data-line-number="924"></td>
        <td id="LC924" class="blob-code js-file-line">            <span class="pl-k">int</span> rcDbz = <span class="pl-c1">mapDataFromRave</span>(scan, &amp;dbzMeta, &amp;dbzImage[<span class="pl-c1">0</span>],<span class="pl-s"><span class="pl-pds">&quot;</span>DBZH<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L925" class="blob-num js-line-number" data-line-number="925"></td>
        <td id="LC925" class="blob-code js-file-line">            <span class="pl-k">if</span> (rcDbz != <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L926" class="blob-num js-line-number" data-line-number="926"></td>
        <td id="LC926" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong while mapping DBZH data from RAVE to LIBVOL2BIRD.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L927" class="blob-num js-line-number" data-line-number="927"></td>
        <td id="LC927" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L928" class="blob-num js-line-number" data-line-number="928"></td>
        <td id="LC928" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L929" class="blob-num js-line-number" data-line-number="929"></td>
        <td id="LC929" class="blob-code js-file-line">            <span class="pl-c">// populate the vradMeta and vradImage variables with data from  </span></td>
      </tr>
      <tr>
        <td id="L930" class="blob-num js-line-number" data-line-number="930"></td>
        <td id="LC930" class="blob-code js-file-line">            <span class="pl-c">// the Rave scan object:</span></td>
      </tr>
      <tr>
        <td id="L931" class="blob-num js-line-number" data-line-number="931"></td>
        <td id="LC931" class="blob-code js-file-line">            <span class="pl-k">int</span> rcVrad = <span class="pl-c1">mapDataFromRave</span>(scan, &amp;vradMeta, &amp;vradImage[<span class="pl-c1">0</span>],<span class="pl-s"><span class="pl-pds">&quot;</span>VRAD<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L932" class="blob-num js-line-number" data-line-number="932"></td>
        <td id="LC932" class="blob-code js-file-line">            <span class="pl-k">if</span> (rcVrad != <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L933" class="blob-num js-line-number" data-line-number="933"></td>
        <td id="LC933" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong while mapping VRAD data from RAVE to LIBVOL2BIRD.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L934" class="blob-num js-line-number" data-line-number="934"></td>
        <td id="LC934" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L935" class="blob-num js-line-number" data-line-number="935"></td>
        <td id="LC935" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L936" class="blob-num js-line-number" data-line-number="936"></td>
        <td id="LC936" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L937" class="blob-num js-line-number" data-line-number="937"></td>
        <td id="LC937" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L938" class="blob-num js-line-number" data-line-number="938"></td>
        <td id="LC938" class="blob-code js-file-line">            <span class="pl-c">//                      calculate vrad texture                   //</span></td>
      </tr>
      <tr>
        <td id="L939" class="blob-num js-line-number" data-line-number="939"></td>
        <td id="LC939" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L940" class="blob-num js-line-number" data-line-number="940"></td>
        <td id="LC940" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L941" class="blob-num js-line-number" data-line-number="941"></td>
        <td id="LC941" class="blob-code js-file-line">            <span class="pl-c1">calcTexture</span>(&amp;texImage[<span class="pl-c1">0</span>], &amp;vradImage[<span class="pl-c1">0</span>], &amp;dbzImage[<span class="pl-c1">0</span>], </td>
      </tr>
      <tr>
        <td id="L942" class="blob-num js-line-number" data-line-number="942"></td>
        <td id="LC942" class="blob-code js-file-line">                        &amp;texMeta, &amp;vradMeta, &amp;dbzMeta);</td>
      </tr>
      <tr>
        <td id="L943" class="blob-num js-line-number" data-line-number="943"></td>
        <td id="LC943" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L944" class="blob-num js-line-number" data-line-number="944"></td>
        <td id="LC944" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L945" class="blob-num js-line-number" data-line-number="945"></td>
        <td id="LC945" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L946" class="blob-num js-line-number" data-line-number="946"></td>
        <td id="LC946" class="blob-code js-file-line">            <span class="pl-c">//        find (weather) cells in the reflectivity image         //</span></td>
      </tr>
      <tr>
        <td id="L947" class="blob-num js-line-number" data-line-number="947"></td>
        <td id="LC947" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L948" class="blob-num js-line-number" data-line-number="948"></td>
        <td id="LC948" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L949" class="blob-num js-line-number" data-line-number="949"></td>
        <td id="LC949" class="blob-code js-file-line">            <span class="pl-k">int</span> nCells = <span class="pl-c1">findWeatherCells</span>(&amp;dbzImage[<span class="pl-c1">0</span>], &amp;cellImage[<span class="pl-c1">0</span>], </td>
      </tr>
      <tr>
        <td id="L950" class="blob-num js-line-number" data-line-number="950"></td>
        <td id="LC950" class="blob-code js-file-line">                                   &amp;dbzMeta);</td>
      </tr>
      <tr>
        <td id="L951" class="blob-num js-line-number" data-line-number="951"></td>
        <td id="LC951" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L952" class="blob-num js-line-number" data-line-number="952"></td>
        <td id="LC952" class="blob-code js-file-line">            <span class="pl-k">if</span> (printCellProp == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L953" class="blob-num js-line-number" data-line-number="953"></td>
        <td id="LC953" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>(<span class="pl-c1">%d</span>/<span class="pl-c1">%d</span>): found <span class="pl-c1">%d</span> cells.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iScan, nScans, nCells);</td>
      </tr>
      <tr>
        <td id="L954" class="blob-num js-line-number" data-line-number="954"></td>
        <td id="LC954" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L955" class="blob-num js-line-number" data-line-number="955"></td>
        <td id="LC955" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L956" class="blob-num js-line-number" data-line-number="956"></td>
        <td id="LC956" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L957" class="blob-num js-line-number" data-line-number="957"></td>
        <td id="LC957" class="blob-code js-file-line">            <span class="pl-c">//                      analyze cells                            //</span></td>
      </tr>
      <tr>
        <td id="L958" class="blob-num js-line-number" data-line-number="958"></td>
        <td id="LC958" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L959" class="blob-num js-line-number" data-line-number="959"></td>
        <td id="LC959" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L960" class="blob-num js-line-number" data-line-number="960"></td>
        <td id="LC960" class="blob-code js-file-line">            <span class="pl-k">int</span> nCellsValid = <span class="pl-c1">analyzeCells</span>(&amp;dbzImage[<span class="pl-c1">0</span>], &amp;vradImage[<span class="pl-c1">0</span>], &amp;texImage[<span class="pl-c1">0</span>], </td>
      </tr>
      <tr>
        <td id="L961" class="blob-num js-line-number" data-line-number="961"></td>
        <td id="LC961" class="blob-code js-file-line">                &amp;clutterImage[<span class="pl-c1">0</span>], &amp;cellImage[<span class="pl-c1">0</span>], &amp;dbzMeta, &amp;vradMeta, &amp;texMeta, </td>
      </tr>
      <tr>
        <td id="L962" class="blob-num js-line-number" data-line-number="962"></td>
        <td id="LC962" class="blob-code js-file-line">                &amp;clutterMeta, nCells, </td>
      </tr>
      <tr>
        <td id="L963" class="blob-num js-line-number" data-line-number="963"></td>
        <td id="LC963" class="blob-code js-file-line">                useStaticClutterData);</td>
      </tr>
      <tr>
        <td id="L964" class="blob-num js-line-number" data-line-number="964"></td>
        <td id="LC964" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L965" class="blob-num js-line-number" data-line-number="965"></td>
        <td id="LC965" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L966" class="blob-num js-line-number" data-line-number="966"></td>
        <td id="LC966" class="blob-code js-file-line">            <span class="pl-c">//                     calculate fringe                          //</span></td>
      </tr>
      <tr>
        <td id="L967" class="blob-num js-line-number" data-line-number="967"></td>
        <td id="LC967" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L968" class="blob-num js-line-number" data-line-number="968"></td>
        <td id="LC968" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L969" class="blob-num js-line-number" data-line-number="969"></td>
        <td id="LC969" class="blob-code js-file-line">            <span class="pl-c1">fringeCells</span>(&amp;cellImage[<span class="pl-c1">0</span>], cellMeta.<span class="pl-smi">nRang</span>, cellMeta.<span class="pl-smi">nAzim</span>, </td>
      </tr>
      <tr>
        <td id="L970" class="blob-num js-line-number" data-line-number="970"></td>
        <td id="LC970" class="blob-code js-file-line">                cellMeta.<span class="pl-smi">azimScale</span>, cellMeta.<span class="pl-smi">rangeScale</span>);</td>
      </tr>
      <tr>
        <td id="L971" class="blob-num js-line-number" data-line-number="971"></td>
        <td id="LC971" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L972" class="blob-num js-line-number" data-line-number="972"></td>
        <td id="LC972" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L973" class="blob-num js-line-number" data-line-number="973"></td>
        <td id="LC973" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L974" class="blob-num js-line-number" data-line-number="974"></td>
        <td id="LC974" class="blob-code js-file-line">            <span class="pl-c">//            print selected outputs to stderr                   //</span></td>
      </tr>
      <tr>
        <td id="L975" class="blob-num js-line-number" data-line-number="975"></td>
        <td id="LC975" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L976" class="blob-num js-line-number" data-line-number="976"></td>
        <td id="LC976" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L977" class="blob-num js-line-number" data-line-number="977"></td>
        <td id="LC977" class="blob-code js-file-line">            <span class="pl-k">if</span> (printDbz == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L978" class="blob-num js-line-number" data-line-number="978"></td>
        <td id="LC978" class="blob-code js-file-line">                <span class="pl-c1">printMeta</span>(&amp;dbzMeta,<span class="pl-s"><span class="pl-pds">&quot;</span>dbzMeta<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L979" class="blob-num js-line-number" data-line-number="979"></td>
        <td id="LC979" class="blob-code js-file-line">                <span class="pl-c1">printImageUChar</span>(&amp;dbzMeta,&amp;dbzImage[<span class="pl-c1">0</span>]);</td>
      </tr>
      <tr>
        <td id="L980" class="blob-num js-line-number" data-line-number="980"></td>
        <td id="LC980" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L981" class="blob-num js-line-number" data-line-number="981"></td>
        <td id="LC981" class="blob-code js-file-line">            <span class="pl-k">if</span> (printVrad == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L982" class="blob-num js-line-number" data-line-number="982"></td>
        <td id="LC982" class="blob-code js-file-line">                <span class="pl-c1">printMeta</span>(&amp;vradMeta,<span class="pl-s"><span class="pl-pds">&quot;</span>vradMeta<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L983" class="blob-num js-line-number" data-line-number="983"></td>
        <td id="LC983" class="blob-code js-file-line">                <span class="pl-c1">printImageUChar</span>(&amp;vradMeta,&amp;vradImage[<span class="pl-c1">0</span>]);</td>
      </tr>
      <tr>
        <td id="L984" class="blob-num js-line-number" data-line-number="984"></td>
        <td id="LC984" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L985" class="blob-num js-line-number" data-line-number="985"></td>
        <td id="LC985" class="blob-code js-file-line">            <span class="pl-k">if</span> (printTex == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L986" class="blob-num js-line-number" data-line-number="986"></td>
        <td id="LC986" class="blob-code js-file-line">                <span class="pl-c1">printMeta</span>(&amp;texMeta,<span class="pl-s"><span class="pl-pds">&quot;</span>texMeta<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L987" class="blob-num js-line-number" data-line-number="987"></td>
        <td id="LC987" class="blob-code js-file-line">                <span class="pl-c1">printImageUChar</span>(&amp;texMeta,&amp;texImage[<span class="pl-c1">0</span>]);</td>
      </tr>
      <tr>
        <td id="L988" class="blob-num js-line-number" data-line-number="988"></td>
        <td id="LC988" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L989" class="blob-num js-line-number" data-line-number="989"></td>
        <td id="LC989" class="blob-code js-file-line">            <span class="pl-k">if</span> (printCell == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L990" class="blob-num js-line-number" data-line-number="990"></td>
        <td id="LC990" class="blob-code js-file-line">                <span class="pl-c1">printMeta</span>(&amp;cellMeta,<span class="pl-s"><span class="pl-pds">&quot;</span>cellMeta<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L991" class="blob-num js-line-number" data-line-number="991"></td>
        <td id="LC991" class="blob-code js-file-line">                <span class="pl-c1">printImageInt</span>(&amp;cellMeta,&amp;cellImage[<span class="pl-c1">0</span>]);</td>
      </tr>
      <tr>
        <td id="L992" class="blob-num js-line-number" data-line-number="992"></td>
        <td id="LC992" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L993" class="blob-num js-line-number" data-line-number="993"></td>
        <td id="LC993" class="blob-code js-file-line">            <span class="pl-k">if</span> (printClut == <span class="pl-c1">TRUE</span>) { </td>
      </tr>
      <tr>
        <td id="L994" class="blob-num js-line-number" data-line-number="994"></td>
        <td id="LC994" class="blob-code js-file-line">                <span class="pl-c1">printMeta</span>(&amp;clutterMeta,<span class="pl-s"><span class="pl-pds">&quot;</span>clutterMeta<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L995" class="blob-num js-line-number" data-line-number="995"></td>
        <td id="LC995" class="blob-code js-file-line">                <span class="pl-c1">printImageUChar</span>(&amp;clutterMeta,&amp;clutterImage[<span class="pl-c1">0</span>]);</td>
      </tr>
      <tr>
        <td id="L996" class="blob-num js-line-number" data-line-number="996"></td>
        <td id="LC996" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L997" class="blob-num js-line-number" data-line-number="997"></td>
        <td id="LC997" class="blob-code js-file-line">                        </td>
      </tr>
      <tr>
        <td id="L998" class="blob-num js-line-number" data-line-number="998"></td>
        <td id="LC998" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L999" class="blob-num js-line-number" data-line-number="999"></td>
        <td id="LC999" class="blob-code js-file-line">            <span class="pl-c">//    fill in the appropriate elements in the points array       //</span></td>
      </tr>
      <tr>
        <td id="L1000" class="blob-num js-line-number" data-line-number="1000"></td>
        <td id="LC1000" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1001" class="blob-num js-line-number" data-line-number="1001"></td>
        <td id="LC1001" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1002" class="blob-num js-line-number" data-line-number="1002"></td>
        <td id="LC1002" class="blob-code js-file-line">            <span class="pl-k">int</span> iLayer;</td>
      </tr>
      <tr>
        <td id="L1003" class="blob-num js-line-number" data-line-number="1003"></td>
        <td id="LC1003" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L1004" class="blob-num js-line-number" data-line-number="1004"></td>
        <td id="LC1004" class="blob-code js-file-line">            <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L1005" class="blob-num js-line-number" data-line-number="1005"></td>
        <td id="LC1005" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L1006" class="blob-num js-line-number" data-line-number="1006"></td>
        <td id="LC1006" class="blob-code js-file-line">                <span class="pl-k">float</span> altitudeMin = iLayer * layerThickness;</td>
      </tr>
      <tr>
        <td id="L1007" class="blob-num js-line-number" data-line-number="1007"></td>
        <td id="LC1007" class="blob-code js-file-line">                <span class="pl-k">float</span> altitudeMax = (iLayer + <span class="pl-c1">1</span>) * layerThickness;</td>
      </tr>
      <tr>
        <td id="L1008" class="blob-num js-line-number" data-line-number="1008"></td>
        <td id="LC1008" class="blob-code js-file-line">                <span class="pl-k">int</span> iRowPoints = indexFrom[iLayer] + nPointsWritten[iLayer];</td>
      </tr>
      <tr>
        <td id="L1009" class="blob-num js-line-number" data-line-number="1009"></td>
        <td id="LC1009" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1010" class="blob-num js-line-number" data-line-number="1010"></td>
        <td id="LC1010" class="blob-code js-file-line">                <span class="pl-k">int</span> n = <span class="pl-c1">getListOfSelectedGates</span>(&amp;vradMeta, &amp;vradImage[<span class="pl-c1">0</span>],</td>
      </tr>
      <tr>
        <td id="L1011" class="blob-num js-line-number" data-line-number="1011"></td>
        <td id="LC1011" class="blob-code js-file-line">                    &amp;dbzMeta, &amp;dbzImage[<span class="pl-c1">0</span>], </td>
      </tr>
      <tr>
        <td id="L1012" class="blob-num js-line-number" data-line-number="1012"></td>
        <td id="LC1012" class="blob-code js-file-line">                    &amp;cellImage[<span class="pl-c1">0</span>], </td>
      </tr>
      <tr>
        <td id="L1013" class="blob-num js-line-number" data-line-number="1013"></td>
        <td id="LC1013" class="blob-code js-file-line">                    altitudeMin, altitudeMax, </td>
      </tr>
      <tr>
        <td id="L1014" class="blob-num js-line-number" data-line-number="1014"></td>
        <td id="LC1014" class="blob-code js-file-line">                    &amp;points[<span class="pl-c1">0</span>], iRowPoints);</td>
      </tr>
      <tr>
        <td id="L1015" class="blob-num js-line-number" data-line-number="1015"></td>
        <td id="LC1015" class="blob-code js-file-line">                    </td>
      </tr>
      <tr>
        <td id="L1016" class="blob-num js-line-number" data-line-number="1016"></td>
        <td id="LC1016" class="blob-code js-file-line">                nPointsWritten[iLayer] += n;</td>
      </tr>
      <tr>
        <td id="L1017" class="blob-num js-line-number" data-line-number="1017"></td>
        <td id="LC1017" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1018" class="blob-num js-line-number" data-line-number="1018"></td>
        <td id="LC1018" class="blob-code js-file-line">                <span class="pl-k">if</span> (indexFrom[iLayer] + nPointsWritten[iLayer] &gt; indexTo[iLayer]) {</td>
      </tr>
      <tr>
        <td id="L1019" class="blob-num js-line-number" data-line-number="1019"></td>
        <td id="LC1019" class="blob-code js-file-line">                    <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Problem occurred: writing over existing data<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1020" class="blob-num js-line-number" data-line-number="1020"></td>
        <td id="LC1020" class="blob-code js-file-line">                    <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L1021" class="blob-num js-line-number" data-line-number="1021"></td>
        <td id="LC1021" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1022" class="blob-num js-line-number" data-line-number="1022"></td>
        <td id="LC1022" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1023" class="blob-num js-line-number" data-line-number="1023"></td>
        <td id="LC1023" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1024" class="blob-num js-line-number" data-line-number="1024"></td>
        <td id="LC1024" class="blob-code js-file-line">            } <span class="pl-c">// endfor (iLayer = 0; iLayer &lt; nLayers; iLayer++)</span></td>
      </tr>
      <tr>
        <td id="L1025" class="blob-num js-line-number" data-line-number="1025"></td>
        <td id="LC1025" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1026" class="blob-num js-line-number" data-line-number="1026"></td>
        <td id="LC1026" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1027" class="blob-num js-line-number" data-line-number="1027"></td>
        <td id="LC1027" class="blob-code js-file-line">            <span class="pl-c">//                         clean up                              //</span></td>
      </tr>
      <tr>
        <td id="L1028" class="blob-num js-line-number" data-line-number="1028"></td>
        <td id="LC1028" class="blob-code js-file-line">            <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1029" class="blob-num js-line-number" data-line-number="1029"></td>
        <td id="LC1029" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1030" class="blob-num js-line-number" data-line-number="1030"></td>
        <td id="LC1030" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1031" class="blob-num js-line-number" data-line-number="1031"></td>
        <td id="LC1031" class="blob-code js-file-line">            <span class="pl-c">// free previously malloc&#39;ed arrays</span></td>
      </tr>
      <tr>
        <td id="L1032" class="blob-num js-line-number" data-line-number="1032"></td>
        <td id="LC1032" class="blob-code js-file-line">            <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) dbzImage);</td>
      </tr>
      <tr>
        <td id="L1033" class="blob-num js-line-number" data-line-number="1033"></td>
        <td id="LC1033" class="blob-code js-file-line">            <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) vradImage);</td>
      </tr>
      <tr>
        <td id="L1034" class="blob-num js-line-number" data-line-number="1034"></td>
        <td id="LC1034" class="blob-code js-file-line">            <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) texImage);</td>
      </tr>
      <tr>
        <td id="L1035" class="blob-num js-line-number" data-line-number="1035"></td>
        <td id="LC1035" class="blob-code js-file-line">            <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) cellImage);</td>
      </tr>
      <tr>
        <td id="L1036" class="blob-num js-line-number" data-line-number="1036"></td>
        <td id="LC1036" class="blob-code js-file-line">            <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) clutterImage);</td>
      </tr>
      <tr>
        <td id="L1037" class="blob-num js-line-number" data-line-number="1037"></td>
        <td id="LC1037" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L1038" class="blob-num js-line-number" data-line-number="1038"></td>
        <td id="LC1038" class="blob-code js-file-line">            <span class="pl-c1">RAVE_OBJECT_RELEASE</span>(scan);</td>
      </tr>
      <tr>
        <td id="L1039" class="blob-num js-line-number" data-line-number="1039"></td>
        <td id="LC1039" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L1040" class="blob-num js-line-number" data-line-number="1040"></td>
        <td id="LC1040" class="blob-code js-file-line">        } <span class="pl-c">// endfor (iScan = 0; iScan &lt; nScans; iScan++)</span></td>
      </tr>
      <tr>
        <td id="L1041" class="blob-num js-line-number" data-line-number="1041"></td>
        <td id="LC1041" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1042" class="blob-num js-line-number" data-line-number="1042"></td>
        <td id="LC1042" class="blob-code js-file-line">}</td>
      </tr>
      <tr>
        <td id="L1043" class="blob-num js-line-number" data-line-number="1043"></td>
        <td id="LC1043" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1044" class="blob-num js-line-number" data-line-number="1044"></td>
        <td id="LC1044" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1045" class="blob-num js-line-number" data-line-number="1045"></td>
        <td id="LC1045" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1046" class="blob-num js-line-number" data-line-number="1046"></td>
        <td id="LC1046" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1047" class="blob-num js-line-number" data-line-number="1047"></td>
        <td id="LC1047" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1048" class="blob-num js-line-number" data-line-number="1048"></td>
        <td id="LC1048" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">detNumberOfGates</span>(<span class="pl-k">const</span> <span class="pl-k">int</span> iLayer, </td>
      </tr>
      <tr>
        <td id="L1049" class="blob-num js-line-number" data-line-number="1049"></td>
        <td id="LC1049" class="blob-code js-file-line">                     <span class="pl-k">const</span> <span class="pl-k">float</span> rangeScale, <span class="pl-k">const</span> <span class="pl-k">float</span> elevAngle,</td>
      </tr>
      <tr>
        <td id="L1050" class="blob-num js-line-number" data-line-number="1050"></td>
        <td id="LC1050" class="blob-code js-file-line">                     <span class="pl-k">const</span> <span class="pl-k">int</span> nRang, <span class="pl-k">const</span> <span class="pl-k">int</span> nAzim,</td>
      </tr>
      <tr>
        <td id="L1051" class="blob-num js-line-number" data-line-number="1051"></td>
        <td id="LC1051" class="blob-code js-file-line">                     <span class="pl-k">const</span> <span class="pl-k">float</span> radarHeight) {</td>
      </tr>
      <tr>
        <td id="L1052" class="blob-num js-line-number" data-line-number="1052"></td>
        <td id="LC1052" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1053" class="blob-num js-line-number" data-line-number="1053"></td>
        <td id="LC1053" class="blob-code js-file-line">    <span class="pl-c">// Determine the number of gates that are within the limits set</span></td>
      </tr>
      <tr>
        <td id="L1054" class="blob-num js-line-number" data-line-number="1054"></td>
        <td id="LC1054" class="blob-code js-file-line">    <span class="pl-c">// by (rangeMin,rangeMax) as well as by (iLayer*layerThickness,</span></td>
      </tr>
      <tr>
        <td id="L1055" class="blob-num js-line-number" data-line-number="1055"></td>
        <td id="LC1055" class="blob-code js-file-line">    <span class="pl-c">// (iLayer+1)*layerThickness).</span></td>
      </tr>
      <tr>
        <td id="L1056" class="blob-num js-line-number" data-line-number="1056"></td>
        <td id="LC1056" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1057" class="blob-num js-line-number" data-line-number="1057"></td>
        <td id="LC1057" class="blob-code js-file-line">    <span class="pl-k">int</span> nGates;</td>
      </tr>
      <tr>
        <td id="L1058" class="blob-num js-line-number" data-line-number="1058"></td>
        <td id="LC1058" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L1059" class="blob-num js-line-number" data-line-number="1059"></td>
        <td id="LC1059" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1060" class="blob-num js-line-number" data-line-number="1060"></td>
        <td id="LC1060" class="blob-code js-file-line">    <span class="pl-k">float</span> layerHeight;</td>
      </tr>
      <tr>
        <td id="L1061" class="blob-num js-line-number" data-line-number="1061"></td>
        <td id="LC1061" class="blob-code js-file-line">    <span class="pl-k">float</span> range;</td>
      </tr>
      <tr>
        <td id="L1062" class="blob-num js-line-number" data-line-number="1062"></td>
        <td id="LC1062" class="blob-code js-file-line">    <span class="pl-k">float</span> beamHeight;</td>
      </tr>
      <tr>
        <td id="L1063" class="blob-num js-line-number" data-line-number="1063"></td>
        <td id="LC1063" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1064" class="blob-num js-line-number" data-line-number="1064"></td>
        <td id="LC1064" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1065" class="blob-num js-line-number" data-line-number="1065"></td>
        <td id="LC1065" class="blob-code js-file-line">    nGates = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1066" class="blob-num js-line-number" data-line-number="1066"></td>
        <td id="LC1066" class="blob-code js-file-line">    layerHeight = (iLayer + <span class="pl-c1">0.5</span>) * layerThickness;</td>
      </tr>
      <tr>
        <td id="L1067" class="blob-num js-line-number" data-line-number="1067"></td>
        <td id="LC1067" class="blob-code js-file-line">    <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L1068" class="blob-num js-line-number" data-line-number="1068"></td>
        <td id="LC1068" class="blob-code js-file-line">        range = (iRang + <span class="pl-c1">0.5</span>) * rangeScale;</td>
      </tr>
      <tr>
        <td id="L1069" class="blob-num js-line-number" data-line-number="1069"></td>
        <td id="LC1069" class="blob-code js-file-line">        <span class="pl-k">if</span> (range &lt; rangeMin || range &gt; rangeMax) {</td>
      </tr>
      <tr>
        <td id="L1070" class="blob-num js-line-number" data-line-number="1070"></td>
        <td id="LC1070" class="blob-code js-file-line">            <span class="pl-c">// the gate is too close to the radar, or too far away</span></td>
      </tr>
      <tr>
        <td id="L1071" class="blob-num js-line-number" data-line-number="1071"></td>
        <td id="LC1071" class="blob-code js-file-line">            <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1072" class="blob-num js-line-number" data-line-number="1072"></td>
        <td id="LC1072" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1073" class="blob-num js-line-number" data-line-number="1073"></td>
        <td id="LC1073" class="blob-code js-file-line">        beamHeight = range * <span class="pl-c1">sin</span>(elevAngle * DEG2RAD) + radarHeight;</td>
      </tr>
      <tr>
        <td id="L1074" class="blob-num js-line-number" data-line-number="1074"></td>
        <td id="LC1074" class="blob-code js-file-line">        <span class="pl-k">if</span> (<span class="pl-c1">fabs</span>(layerHeight - beamHeight) &gt; <span class="pl-c1">0.5</span>*layerThickness) {</td>
      </tr>
      <tr>
        <td id="L1075" class="blob-num js-line-number" data-line-number="1075"></td>
        <td id="LC1075" class="blob-code js-file-line">            <span class="pl-c">// the gate is not close enough to the altitude layer of interest</span></td>
      </tr>
      <tr>
        <td id="L1076" class="blob-num js-line-number" data-line-number="1076"></td>
        <td id="LC1076" class="blob-code js-file-line">            <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1077" class="blob-num js-line-number" data-line-number="1077"></td>
        <td id="LC1077" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1078" class="blob-num js-line-number" data-line-number="1078"></td>
        <td id="LC1078" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1079" class="blob-num js-line-number" data-line-number="1079"></td>
        <td id="LC1079" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1080" class="blob-num js-line-number" data-line-number="1080"></td>
        <td id="LC1080" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>iRang = <span class="pl-c1">%d</span>; range = <span class="pl-c1">%f</span>; beamHeight = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iRang,range,beamHeight);</td>
      </tr>
      <tr>
        <td id="L1081" class="blob-num js-line-number" data-line-number="1081"></td>
        <td id="LC1081" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1082" class="blob-num js-line-number" data-line-number="1082"></td>
        <td id="LC1082" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1083" class="blob-num js-line-number" data-line-number="1083"></td>
        <td id="LC1083" class="blob-code js-file-line">        nGates += nAzim;</td>
      </tr>
      <tr>
        <td id="L1084" class="blob-num js-line-number" data-line-number="1084"></td>
        <td id="LC1084" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1085" class="blob-num js-line-number" data-line-number="1085"></td>
        <td id="LC1085" class="blob-code js-file-line">    } <span class="pl-c">// for iRang</span></td>
      </tr>
      <tr>
        <td id="L1086" class="blob-num js-line-number" data-line-number="1086"></td>
        <td id="LC1086" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1087" class="blob-num js-line-number" data-line-number="1087"></td>
        <td id="LC1087" class="blob-code js-file-line">    <span class="pl-k">return</span> nGates;</td>
      </tr>
      <tr>
        <td id="L1088" class="blob-num js-line-number" data-line-number="1088"></td>
        <td id="LC1088" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1089" class="blob-num js-line-number" data-line-number="1089"></td>
        <td id="LC1089" class="blob-code js-file-line">} <span class="pl-c">// detNumberOfGates()</span></td>
      </tr>
      <tr>
        <td id="L1090" class="blob-num js-line-number" data-line-number="1090"></td>
        <td id="LC1090" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1091" class="blob-num js-line-number" data-line-number="1091"></td>
        <td id="LC1091" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1092" class="blob-num js-line-number" data-line-number="1092"></td>
        <td id="LC1092" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1093" class="blob-num js-line-number" data-line-number="1093"></td>
        <td id="LC1093" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1094" class="blob-num js-line-number" data-line-number="1094"></td>
        <td id="LC1094" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1095" class="blob-num js-line-number" data-line-number="1095"></td>
        <td id="LC1095" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">detSvdfitArraySize</span>(PolarVolume_t* volume) {</td>
      </tr>
      <tr>
        <td id="L1096" class="blob-num js-line-number" data-line-number="1096"></td>
        <td id="LC1096" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1097" class="blob-num js-line-number" data-line-number="1097"></td>
        <td id="LC1097" class="blob-code js-file-line">    <span class="pl-k">int</span> iScan;</td>
      </tr>
      <tr>
        <td id="L1098" class="blob-num js-line-number" data-line-number="1098"></td>
        <td id="LC1098" class="blob-code js-file-line">    <span class="pl-k">int</span> nScans = <span class="pl-c1">PolarVolume_getNumberOfScans</span>(volume);</td>
      </tr>
      <tr>
        <td id="L1099" class="blob-num js-line-number" data-line-number="1099"></td>
        <td id="LC1099" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1100" class="blob-num js-line-number" data-line-number="1100"></td>
        <td id="LC1100" class="blob-code js-file-line">    <span class="pl-k">int</span> iLayer;</td>
      </tr>
      <tr>
        <td id="L1101" class="blob-num js-line-number" data-line-number="1101"></td>
        <td id="LC1101" class="blob-code js-file-line">    <span class="pl-k">int</span> nRowsPoints = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1102" class="blob-num js-line-number" data-line-number="1102"></td>
        <td id="LC1102" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1103" class="blob-num js-line-number" data-line-number="1103"></td>
        <td id="LC1103" class="blob-code js-file-line">    <span class="pl-k">int</span>* nGates = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nLayers);</td>
      </tr>
      <tr>
        <td id="L1104" class="blob-num js-line-number" data-line-number="1104"></td>
        <td id="LC1104" class="blob-code js-file-line">    <span class="pl-k">int</span>* nGatesAcc = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nLayers);</td>
      </tr>
      <tr>
        <td id="L1105" class="blob-num js-line-number" data-line-number="1105"></td>
        <td id="LC1105" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1106" class="blob-num js-line-number" data-line-number="1106"></td>
        <td id="LC1106" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L1107" class="blob-num js-line-number" data-line-number="1107"></td>
        <td id="LC1107" class="blob-code js-file-line">        nGates[iLayer] = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1108" class="blob-num js-line-number" data-line-number="1108"></td>
        <td id="LC1108" class="blob-code js-file-line">        nGatesAcc[iLayer] = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1109" class="blob-num js-line-number" data-line-number="1109"></td>
        <td id="LC1109" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1110" class="blob-num js-line-number" data-line-number="1110"></td>
        <td id="LC1110" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1111" class="blob-num js-line-number" data-line-number="1111"></td>
        <td id="LC1111" class="blob-code js-file-line">    <span class="pl-k">for</span> (iScan = <span class="pl-c1">0</span>; iScan &lt; nScans; iScan++) {</td>
      </tr>
      <tr>
        <td id="L1112" class="blob-num js-line-number" data-line-number="1112"></td>
        <td id="LC1112" class="blob-code js-file-line">        <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L1113" class="blob-num js-line-number" data-line-number="1113"></td>
        <td id="LC1113" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1114" class="blob-num js-line-number" data-line-number="1114"></td>
        <td id="LC1114" class="blob-code js-file-line">            PolarScan_t* scan = <span class="pl-c1">PolarVolume_getScan</span>(volume, iScan);</td>
      </tr>
      <tr>
        <td id="L1115" class="blob-num js-line-number" data-line-number="1115"></td>
        <td id="LC1115" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L1116" class="blob-num js-line-number" data-line-number="1116"></td>
        <td id="LC1116" class="blob-code js-file-line">            <span class="pl-k">int</span> nRang = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNbins</span>(scan);</td>
      </tr>
      <tr>
        <td id="L1117" class="blob-num js-line-number" data-line-number="1117"></td>
        <td id="LC1117" class="blob-code js-file-line">            <span class="pl-k">int</span> nAzim = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNrays</span>(scan);</td>
      </tr>
      <tr>
        <td id="L1118" class="blob-num js-line-number" data-line-number="1118"></td>
        <td id="LC1118" class="blob-code js-file-line">            <span class="pl-k">float</span> elevAngle = (<span class="pl-k">float</span>) (<span class="pl-c1">360</span> * <span class="pl-c1">PolarScan_getElangle</span>(scan) / <span class="pl-c1">2</span> / PI );</td>
      </tr>
      <tr>
        <td id="L1119" class="blob-num js-line-number" data-line-number="1119"></td>
        <td id="LC1119" class="blob-code js-file-line">            <span class="pl-k">float</span> rangeScale = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScan_getRscale</span>(scan);</td>
      </tr>
      <tr>
        <td id="L1120" class="blob-num js-line-number" data-line-number="1120"></td>
        <td id="LC1120" class="blob-code js-file-line">            <span class="pl-k">float</span> radarHeight = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScan_getHeight</span>(scan);</td>
      </tr>
      <tr>
        <td id="L1121" class="blob-num js-line-number" data-line-number="1121"></td>
        <td id="LC1121" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1122" class="blob-num js-line-number" data-line-number="1122"></td>
        <td id="LC1122" class="blob-code js-file-line">            nGates[iLayer] += <span class="pl-c1">detNumberOfGates</span>(iLayer, </td>
      </tr>
      <tr>
        <td id="L1123" class="blob-num js-line-number" data-line-number="1123"></td>
        <td id="LC1123" class="blob-code js-file-line">                rangeScale, elevAngle, nRang, </td>
      </tr>
      <tr>
        <td id="L1124" class="blob-num js-line-number" data-line-number="1124"></td>
        <td id="LC1124" class="blob-code js-file-line">                nAzim, radarHeight);</td>
      </tr>
      <tr>
        <td id="L1125" class="blob-num js-line-number" data-line-number="1125"></td>
        <td id="LC1125" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L1126" class="blob-num js-line-number" data-line-number="1126"></td>
        <td id="LC1126" class="blob-code js-file-line">            <span class="pl-c1">RAVE_OBJECT_RELEASE</span>(scan);</td>
      </tr>
      <tr>
        <td id="L1127" class="blob-num js-line-number" data-line-number="1127"></td>
        <td id="LC1127" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1128" class="blob-num js-line-number" data-line-number="1128"></td>
        <td id="LC1128" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1129" class="blob-num js-line-number" data-line-number="1129"></td>
        <td id="LC1129" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1130" class="blob-num js-line-number" data-line-number="1130"></td>
        <td id="LC1130" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1131" class="blob-num js-line-number" data-line-number="1131"></td>
        <td id="LC1131" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L1132" class="blob-num js-line-number" data-line-number="1132"></td>
        <td id="LC1132" class="blob-code js-file-line">        <span class="pl-k">if</span> (iLayer == <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L1133" class="blob-num js-line-number" data-line-number="1133"></td>
        <td id="LC1133" class="blob-code js-file-line">            nGatesAcc[iLayer] = nGates[iLayer];</td>
      </tr>
      <tr>
        <td id="L1134" class="blob-num js-line-number" data-line-number="1134"></td>
        <td id="LC1134" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1135" class="blob-num js-line-number" data-line-number="1135"></td>
        <td id="LC1135" class="blob-code js-file-line">        <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L1136" class="blob-num js-line-number" data-line-number="1136"></td>
        <td id="LC1136" class="blob-code js-file-line">            nGatesAcc[iLayer] = nGatesAcc[iLayer-<span class="pl-c1">1</span>] + nGates[iLayer];</td>
      </tr>
      <tr>
        <td id="L1137" class="blob-num js-line-number" data-line-number="1137"></td>
        <td id="LC1137" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1138" class="blob-num js-line-number" data-line-number="1138"></td>
        <td id="LC1138" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1139" class="blob-num js-line-number" data-line-number="1139"></td>
        <td id="LC1139" class="blob-code js-file-line">    nRowsPoints = nGatesAcc[nLayers-<span class="pl-c1">1</span>];</td>
      </tr>
      <tr>
        <td id="L1140" class="blob-num js-line-number" data-line-number="1140"></td>
        <td id="LC1140" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1141" class="blob-num js-line-number" data-line-number="1141"></td>
        <td id="LC1141" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L1142" class="blob-num js-line-number" data-line-number="1142"></td>
        <td id="LC1142" class="blob-code js-file-line">        <span class="pl-k">if</span> (iLayer == <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L1143" class="blob-num js-line-number" data-line-number="1143"></td>
        <td id="LC1143" class="blob-code js-file-line">            indexFrom[iLayer] = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1144" class="blob-num js-line-number" data-line-number="1144"></td>
        <td id="LC1144" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1145" class="blob-num js-line-number" data-line-number="1145"></td>
        <td id="LC1145" class="blob-code js-file-line">        <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L1146" class="blob-num js-line-number" data-line-number="1146"></td>
        <td id="LC1146" class="blob-code js-file-line">            indexFrom[iLayer] = nGatesAcc[iLayer-<span class="pl-c1">1</span>];</td>
      </tr>
      <tr>
        <td id="L1147" class="blob-num js-line-number" data-line-number="1147"></td>
        <td id="LC1147" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1148" class="blob-num js-line-number" data-line-number="1148"></td>
        <td id="LC1148" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1149" class="blob-num js-line-number" data-line-number="1149"></td>
        <td id="LC1149" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L1150" class="blob-num js-line-number" data-line-number="1150"></td>
        <td id="LC1150" class="blob-code js-file-line">        indexTo[iLayer] = nGatesAcc[iLayer];</td>
      </tr>
      <tr>
        <td id="L1151" class="blob-num js-line-number" data-line-number="1151"></td>
        <td id="LC1151" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1152" class="blob-num js-line-number" data-line-number="1152"></td>
        <td id="LC1152" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1153" class="blob-num js-line-number" data-line-number="1153"></td>
        <td id="LC1153" class="blob-code js-file-line">    <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) nGates);</td>
      </tr>
      <tr>
        <td id="L1154" class="blob-num js-line-number" data-line-number="1154"></td>
        <td id="LC1154" class="blob-code js-file-line">    <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) nGatesAcc);</td>
      </tr>
      <tr>
        <td id="L1155" class="blob-num js-line-number" data-line-number="1155"></td>
        <td id="LC1155" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1156" class="blob-num js-line-number" data-line-number="1156"></td>
        <td id="LC1156" class="blob-code js-file-line">    <span class="pl-k">return</span> nRowsPoints;</td>
      </tr>
      <tr>
        <td id="L1157" class="blob-num js-line-number" data-line-number="1157"></td>
        <td id="LC1157" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1158" class="blob-num js-line-number" data-line-number="1158"></td>
        <td id="LC1158" class="blob-code js-file-line">}  <span class="pl-c">// detSvdfitArraySize()</span></td>
      </tr>
      <tr>
        <td id="L1159" class="blob-num js-line-number" data-line-number="1159"></td>
        <td id="LC1159" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1160" class="blob-num js-line-number" data-line-number="1160"></td>
        <td id="LC1160" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1161" class="blob-num js-line-number" data-line-number="1161"></td>
        <td id="LC1161" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1162" class="blob-num js-line-number" data-line-number="1162"></td>
        <td id="LC1162" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1163" class="blob-num js-line-number" data-line-number="1163"></td>
        <td id="LC1163" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1164" class="blob-num js-line-number" data-line-number="1164"></td>
        <td id="LC1164" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1165" class="blob-num js-line-number" data-line-number="1165"></td>
        <td id="LC1165" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">findWeatherCells</span>(<span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *dbzImage, <span class="pl-k">int</span> *cellImage,</td>
      </tr>
      <tr>
        <td id="L1166" class="blob-num js-line-number" data-line-number="1166"></td>
        <td id="LC1166" class="blob-code js-file-line">        <span class="pl-k">const</span> SCANMETA *dbzMeta) {</td>
      </tr>
      <tr>
        <td id="L1167" class="blob-num js-line-number" data-line-number="1167"></td>
        <td id="LC1167" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1168" class="blob-num js-line-number" data-line-number="1168"></td>
        <td id="LC1168" class="blob-code js-file-line">    <span class="pl-c">//  ----------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1169" class="blob-num js-line-number" data-line-number="1169"></td>
        <td id="LC1169" class="blob-code js-file-line">    <span class="pl-c">//  This function detects the cells in &#39;dbzImage&#39; using a threshold value of      //</span></td>
      </tr>
      <tr>
        <td id="L1170" class="blob-num js-line-number" data-line-number="1170"></td>
        <td id="LC1170" class="blob-code js-file-line">    <span class="pl-c">//  &#39;dbzThresMin&#39; and a non-recursive algorithm which looks for neighboring       //</span></td>
      </tr>
      <tr>
        <td id="L1171" class="blob-num js-line-number" data-line-number="1171"></td>
        <td id="LC1171" class="blob-code js-file-line">    <span class="pl-c">//  pixels above that threshold. On return the marked cells are contained by      //</span></td>
      </tr>
      <tr>
        <td id="L1172" class="blob-num js-line-number" data-line-number="1172"></td>
        <td id="LC1172" class="blob-code js-file-line">    <span class="pl-c">//  &#39;cellImage&#39;. The number of detected cells/highest index value is returned.    //</span></td>
      </tr>
      <tr>
        <td id="L1173" class="blob-num js-line-number" data-line-number="1173"></td>
        <td id="LC1173" class="blob-code js-file-line">    <span class="pl-c">//  ----------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1174" class="blob-num js-line-number" data-line-number="1174"></td>
        <td id="LC1174" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1175" class="blob-num js-line-number" data-line-number="1175"></td>
        <td id="LC1175" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1176" class="blob-num js-line-number" data-line-number="1176"></td>
        <td id="LC1176" class="blob-code js-file-line">    <span class="pl-k">int</span> iCellIdentifier;</td>
      </tr>
      <tr>
        <td id="L1177" class="blob-num js-line-number" data-line-number="1177"></td>
        <td id="LC1177" class="blob-code js-file-line">    <span class="pl-k">int</span> nCells;</td>
      </tr>
      <tr>
        <td id="L1178" class="blob-num js-line-number" data-line-number="1178"></td>
        <td id="LC1178" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L1179" class="blob-num js-line-number" data-line-number="1179"></td>
        <td id="LC1179" class="blob-code js-file-line">    <span class="pl-k">int</span> nRang;</td>
      </tr>
      <tr>
        <td id="L1180" class="blob-num js-line-number" data-line-number="1180"></td>
        <td id="LC1180" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzim;</td>
      </tr>
      <tr>
        <td id="L1181" class="blob-num js-line-number" data-line-number="1181"></td>
        <td id="LC1181" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzim;</td>
      </tr>
      <tr>
        <td id="L1182" class="blob-num js-line-number" data-line-number="1182"></td>
        <td id="LC1182" class="blob-code js-file-line">    <span class="pl-k">int</span> iNeighborhood;</td>
      </tr>
      <tr>
        <td id="L1183" class="blob-num js-line-number" data-line-number="1183"></td>
        <td id="LC1183" class="blob-code js-file-line">    <span class="pl-k">int</span> nNeighborhood;</td>
      </tr>
      <tr>
        <td id="L1184" class="blob-num js-line-number" data-line-number="1184"></td>
        <td id="LC1184" class="blob-code js-file-line">    <span class="pl-k">int</span> count;</td>
      </tr>
      <tr>
        <td id="L1185" class="blob-num js-line-number" data-line-number="1185"></td>
        <td id="LC1185" class="blob-code js-file-line">    <span class="pl-k">int</span> cellImageInitialValue;</td>
      </tr>
      <tr>
        <td id="L1186" class="blob-num js-line-number" data-line-number="1186"></td>
        <td id="LC1186" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1187" class="blob-num js-line-number" data-line-number="1187"></td>
        <td id="LC1187" class="blob-code js-file-line">    <span class="pl-k">unsigned</span> <span class="pl-k">char</span> dbzThres;</td>
      </tr>
      <tr>
        <td id="L1188" class="blob-num js-line-number" data-line-number="1188"></td>
        <td id="LC1188" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1189" class="blob-num js-line-number" data-line-number="1189"></td>
        <td id="LC1189" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L1190" class="blob-num js-line-number" data-line-number="1190"></td>
        <td id="LC1190" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobalOther;</td>
      </tr>
      <tr>
        <td id="L1191" class="blob-num js-line-number" data-line-number="1191"></td>
        <td id="LC1191" class="blob-code js-file-line">    <span class="pl-k">int</span> nGlobal;</td>
      </tr>
      <tr>
        <td id="L1192" class="blob-num js-line-number" data-line-number="1192"></td>
        <td id="LC1192" class="blob-code js-file-line">    <span class="pl-k">int</span> iLocal;</td>
      </tr>
      <tr>
        <td id="L1193" class="blob-num js-line-number" data-line-number="1193"></td>
        <td id="LC1193" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1194" class="blob-num js-line-number" data-line-number="1194"></td>
        <td id="LC1194" class="blob-code js-file-line">    <span class="pl-k">unsigned</span> <span class="pl-k">char</span> dbzMissing;</td>
      </tr>
      <tr>
        <td id="L1195" class="blob-num js-line-number" data-line-number="1195"></td>
        <td id="LC1195" class="blob-code js-file-line">    <span class="pl-k">int</span> dbznAzim;</td>
      </tr>
      <tr>
        <td id="L1196" class="blob-num js-line-number" data-line-number="1196"></td>
        <td id="LC1196" class="blob-code js-file-line">    <span class="pl-k">int</span> dbznRang;</td>
      </tr>
      <tr>
        <td id="L1197" class="blob-num js-line-number" data-line-number="1197"></td>
        <td id="LC1197" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzValueOffset;</td>
      </tr>
      <tr>
        <td id="L1198" class="blob-num js-line-number" data-line-number="1198"></td>
        <td id="LC1198" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzValueScale;</td>
      </tr>
      <tr>
        <td id="L1199" class="blob-num js-line-number" data-line-number="1199"></td>
        <td id="LC1199" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzRangeScale;</td>
      </tr>
      <tr>
        <td id="L1200" class="blob-num js-line-number" data-line-number="1200"></td>
        <td id="LC1200" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1201" class="blob-num js-line-number" data-line-number="1201"></td>
        <td id="LC1201" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1202" class="blob-num js-line-number" data-line-number="1202"></td>
        <td id="LC1202" class="blob-code js-file-line">    <span class="pl-c">// These next two local variables &#39;nAzimNeighborhood&#39; and </span></td>
      </tr>
      <tr>
        <td id="L1203" class="blob-num js-line-number" data-line-number="1203"></td>
        <td id="LC1203" class="blob-code js-file-line">    <span class="pl-c">// &#39;nRangNeighborhood&#39; shadow two static global variables of the </span></td>
      </tr>
      <tr>
        <td id="L1204" class="blob-num js-line-number" data-line-number="1204"></td>
        <td id="LC1204" class="blob-code js-file-line">    <span class="pl-c">// same name, but that is in fact what you want. This is because the</span></td>
      </tr>
      <tr>
        <td id="L1205" class="blob-num js-line-number" data-line-number="1205"></td>
        <td id="LC1205" class="blob-code js-file-line">    <span class="pl-c">// value of &#39;nAzimNeighborHood&#39; (function scope) is not equal to </span></td>
      </tr>
      <tr>
        <td id="L1206" class="blob-num js-line-number" data-line-number="1206"></td>
        <td id="LC1206" class="blob-code js-file-line">    <span class="pl-c">// &#39;nAzimNeighborHood&#39; (file scope) and because &#39;nRangNeighborHood&#39; </span></td>
      </tr>
      <tr>
        <td id="L1207" class="blob-num js-line-number" data-line-number="1207"></td>
        <td id="LC1207" class="blob-code js-file-line">    <span class="pl-c">// (function scope) is not equal to &#39;nRangNeighborHood&#39; (file </span></td>
      </tr>
      <tr>
        <td id="L1208" class="blob-num js-line-number" data-line-number="1208"></td>
        <td id="LC1208" class="blob-code js-file-line">    <span class="pl-c">// scope). The variable names are the same though, because they play</span></td>
      </tr>
      <tr>
        <td id="L1209" class="blob-num js-line-number" data-line-number="1209"></td>
        <td id="LC1209" class="blob-code js-file-line">    <span class="pl-c">// the same role in &#39;findNearbyGateIndex()&#39;.</span></td>
      </tr>
      <tr>
        <td id="L1210" class="blob-num js-line-number" data-line-number="1210"></td>
        <td id="LC1210" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzimNeighborhood;</td>
      </tr>
      <tr>
        <td id="L1211" class="blob-num js-line-number" data-line-number="1211"></td>
        <td id="LC1211" class="blob-code js-file-line">    <span class="pl-k">int</span> nRangNeighborhood;</td>
      </tr>
      <tr>
        <td id="L1212" class="blob-num js-line-number" data-line-number="1212"></td>
        <td id="LC1212" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1213" class="blob-num js-line-number" data-line-number="1213"></td>
        <td id="LC1213" class="blob-code js-file-line">    <span class="pl-k">int</span> nHalfNeighborhood;</td>
      </tr>
      <tr>
        <td id="L1214" class="blob-num js-line-number" data-line-number="1214"></td>
        <td id="LC1214" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1215" class="blob-num js-line-number" data-line-number="1215"></td>
        <td id="LC1215" class="blob-code js-file-line">    <span class="pl-k">int</span> cellIdentifierGlobal;</td>
      </tr>
      <tr>
        <td id="L1216" class="blob-num js-line-number" data-line-number="1216"></td>
        <td id="LC1216" class="blob-code js-file-line">    <span class="pl-k">int</span> cellIdentifierGlobalOther;</td>
      </tr>
      <tr>
        <td id="L1217" class="blob-num js-line-number" data-line-number="1217"></td>
        <td id="LC1217" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobalInner;</td>
      </tr>
      <tr>
        <td id="L1218" class="blob-num js-line-number" data-line-number="1218"></td>
        <td id="LC1218" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1219" class="blob-num js-line-number" data-line-number="1219"></td>
        <td id="LC1219" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1220" class="blob-num js-line-number" data-line-number="1220"></td>
        <td id="LC1220" class="blob-code js-file-line">    #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1221" class="blob-num js-line-number" data-line-number="1221"></td>
        <td id="LC1221" class="blob-code js-file-line">    <span class="pl-k">int</span> dbg = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1222" class="blob-num js-line-number" data-line-number="1222"></td>
        <td id="LC1222" class="blob-code js-file-line">    #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1223" class="blob-num js-line-number" data-line-number="1223"></td>
        <td id="LC1223" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1224" class="blob-num js-line-number" data-line-number="1224"></td>
        <td id="LC1224" class="blob-code js-file-line">    <span class="pl-k">if</span> (dbzImage==<span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L1225" class="blob-num js-line-number" data-line-number="1225"></td>
        <td id="LC1225" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Input argument dbzImage is NULL.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1226" class="blob-num js-line-number" data-line-number="1226"></td>
        <td id="LC1226" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1227" class="blob-num js-line-number" data-line-number="1227"></td>
        <td id="LC1227" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1228" class="blob-num js-line-number" data-line-number="1228"></td>
        <td id="LC1228" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1229" class="blob-num js-line-number" data-line-number="1229"></td>
        <td id="LC1229" class="blob-code js-file-line">    dbzMissing = dbzMeta-&gt;missing;</td>
      </tr>
      <tr>
        <td id="L1230" class="blob-num js-line-number" data-line-number="1230"></td>
        <td id="LC1230" class="blob-code js-file-line">    dbznAzim = dbzMeta-&gt;nAzim;</td>
      </tr>
      <tr>
        <td id="L1231" class="blob-num js-line-number" data-line-number="1231"></td>
        <td id="LC1231" class="blob-code js-file-line">    dbznRang = dbzMeta-&gt;nRang;</td>
      </tr>
      <tr>
        <td id="L1232" class="blob-num js-line-number" data-line-number="1232"></td>
        <td id="LC1232" class="blob-code js-file-line">    dbzValueOffset = dbzMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L1233" class="blob-num js-line-number" data-line-number="1233"></td>
        <td id="LC1233" class="blob-code js-file-line">    dbzValueScale = dbzMeta-&gt;valueScale;</td>
      </tr>
      <tr>
        <td id="L1234" class="blob-num js-line-number" data-line-number="1234"></td>
        <td id="LC1234" class="blob-code js-file-line">    dbzRangeScale = dbzMeta-&gt;rangeScale;</td>
      </tr>
      <tr>
        <td id="L1235" class="blob-num js-line-number" data-line-number="1235"></td>
        <td id="LC1235" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1236" class="blob-num js-line-number" data-line-number="1236"></td>
        <td id="LC1236" class="blob-code js-file-line">    nAzim = dbznAzim;</td>
      </tr>
      <tr>
        <td id="L1237" class="blob-num js-line-number" data-line-number="1237"></td>
        <td id="LC1237" class="blob-code js-file-line">    nRang = dbznRang;</td>
      </tr>
      <tr>
        <td id="L1238" class="blob-num js-line-number" data-line-number="1238"></td>
        <td id="LC1238" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1239" class="blob-num js-line-number" data-line-number="1239"></td>
        <td id="LC1239" class="blob-code js-file-line">    nGlobal = nAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L1240" class="blob-num js-line-number" data-line-number="1240"></td>
        <td id="LC1240" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1241" class="blob-num js-line-number" data-line-number="1241"></td>
        <td id="LC1241" class="blob-code js-file-line">    <span class="pl-c">// We use a neighborhood of 3x3, because in this function we are</span></td>
      </tr>
      <tr>
        <td id="L1242" class="blob-num js-line-number" data-line-number="1242"></td>
        <td id="LC1242" class="blob-code js-file-line">    <span class="pl-c">// interested in a cell&#39;s direct neighbors. See also comment at</span></td>
      </tr>
      <tr>
        <td id="L1243" class="blob-num js-line-number" data-line-number="1243"></td>
        <td id="LC1243" class="blob-code js-file-line">    <span class="pl-c">// variable definitions of &#39;nAzimNeighborhood&#39; and &#39;nRangNeighborhood&#39;</span></td>
      </tr>
      <tr>
        <td id="L1244" class="blob-num js-line-number" data-line-number="1244"></td>
        <td id="LC1244" class="blob-code js-file-line">    <span class="pl-c">// above. </span></td>
      </tr>
      <tr>
        <td id="L1245" class="blob-num js-line-number" data-line-number="1245"></td>
        <td id="LC1245" class="blob-code js-file-line">    nAzimNeighborhood = <span class="pl-c1">3</span>;</td>
      </tr>
      <tr>
        <td id="L1246" class="blob-num js-line-number" data-line-number="1246"></td>
        <td id="LC1246" class="blob-code js-file-line">    nRangNeighborhood = <span class="pl-c1">3</span>;</td>
      </tr>
      <tr>
        <td id="L1247" class="blob-num js-line-number" data-line-number="1247"></td>
        <td id="LC1247" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1248" class="blob-num js-line-number" data-line-number="1248"></td>
        <td id="LC1248" class="blob-code js-file-line">    nNeighborhood = nAzimNeighborhood * nRangNeighborhood;</td>
      </tr>
      <tr>
        <td id="L1249" class="blob-num js-line-number" data-line-number="1249"></td>
        <td id="LC1249" class="blob-code js-file-line">    nHalfNeighborhood = (nNeighborhood - <span class="pl-c1">1</span>)/<span class="pl-c1">2</span>;</td>
      </tr>
      <tr>
        <td id="L1250" class="blob-num js-line-number" data-line-number="1250"></td>
        <td id="LC1250" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1251" class="blob-num js-line-number" data-line-number="1251"></td>
        <td id="LC1251" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1252" class="blob-num js-line-number" data-line-number="1252"></td>
        <td id="LC1252" class="blob-code js-file-line">    <span class="pl-k">if</span> (dbzImage != <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L1253" class="blob-num js-line-number" data-line-number="1253"></td>
        <td id="LC1253" class="blob-code js-file-line">        dbzThres = (<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) <span class="pl-c1">ROUND</span>((dbzThresMin - dbzValueOffset) / dbzValueScale);</td>
      </tr>
      <tr>
        <td id="L1254" class="blob-num js-line-number" data-line-number="1254"></td>
        <td id="LC1254" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1255" class="blob-num js-line-number" data-line-number="1255"></td>
        <td id="LC1255" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1256" class="blob-num js-line-number" data-line-number="1256"></td>
        <td id="LC1256" class="blob-code js-file-line">    cellImageInitialValue = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1257" class="blob-num js-line-number" data-line-number="1257"></td>
        <td id="LC1257" class="blob-code js-file-line">    <span class="pl-k">for</span> (iGlobal = <span class="pl-c1">0</span>; iGlobal &lt; nGlobal; iGlobal++) {</td>
      </tr>
      <tr>
        <td id="L1258" class="blob-num js-line-number" data-line-number="1258"></td>
        <td id="LC1258" class="blob-code js-file-line">        cellImage[iGlobal] = cellImageInitialValue;</td>
      </tr>
      <tr>
        <td id="L1259" class="blob-num js-line-number" data-line-number="1259"></td>
        <td id="LC1259" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1260" class="blob-num js-line-number" data-line-number="1260"></td>
        <td id="LC1260" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1261" class="blob-num js-line-number" data-line-number="1261"></td>
        <td id="LC1261" class="blob-code js-file-line">    <span class="pl-c">// If threshold value is equal to missing value, return.</span></td>
      </tr>
      <tr>
        <td id="L1262" class="blob-num js-line-number" data-line-number="1262"></td>
        <td id="LC1262" class="blob-code js-file-line">    <span class="pl-k">if</span> (dbzThres == dbzMissing) {</td>
      </tr>
      <tr>
        <td id="L1263" class="blob-num js-line-number" data-line-number="1263"></td>
        <td id="LC1263" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1264" class="blob-num js-line-number" data-line-number="1264"></td>
        <td id="LC1264" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1265" class="blob-num js-line-number" data-line-number="1265"></td>
        <td id="LC1265" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1266" class="blob-num js-line-number" data-line-number="1266"></td>
        <td id="LC1266" class="blob-code js-file-line">    <span class="pl-c">// ----------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1267" class="blob-num js-line-number" data-line-number="1267"></td>
        <td id="LC1267" class="blob-code js-file-line">    <span class="pl-c">// Labeling of groups of connected pixels using horizontal, vertical, and  //</span></td>
      </tr>
      <tr>
        <td id="L1268" class="blob-num js-line-number" data-line-number="1268"></td>
        <td id="LC1268" class="blob-code js-file-line">    <span class="pl-c">// diagonal connections. The algorithm is described in &#39;Digital Image      //</span></td>
      </tr>
      <tr>
        <td id="L1269" class="blob-num js-line-number" data-line-number="1269"></td>
        <td id="LC1269" class="blob-code js-file-line">    <span class="pl-c">// Processing&#39; by Gonzales and Woods published by Addison-Wesley.          //</span></td>
      </tr>
      <tr>
        <td id="L1270" class="blob-num js-line-number" data-line-number="1270"></td>
        <td id="LC1270" class="blob-code js-file-line">    <span class="pl-c">// ----------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1271" class="blob-num js-line-number" data-line-number="1271"></td>
        <td id="LC1271" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1272" class="blob-num js-line-number" data-line-number="1272"></td>
        <td id="LC1272" class="blob-code js-file-line">    iCellIdentifier = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1273" class="blob-num js-line-number" data-line-number="1273"></td>
        <td id="LC1273" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1274" class="blob-num js-line-number" data-line-number="1274"></td>
        <td id="LC1274" class="blob-code js-file-line">    <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) {</td>
      </tr>
      <tr>
        <td id="L1275" class="blob-num js-line-number" data-line-number="1275"></td>
        <td id="LC1275" class="blob-code js-file-line">        <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L1276" class="blob-num js-line-number" data-line-number="1276"></td>
        <td id="LC1276" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1277" class="blob-num js-line-number" data-line-number="1277"></td>
        <td id="LC1277" class="blob-code js-file-line">            iGlobal = iRang + iAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L1278" class="blob-num js-line-number" data-line-number="1278"></td>
        <td id="LC1278" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1279" class="blob-num js-line-number" data-line-number="1279"></td>
        <td id="LC1279" class="blob-code js-file-line">            <span class="pl-k">if</span> ((<span class="pl-k">float</span>)(iRang + <span class="pl-c1">1</span>) * dbzRangeScale &gt; rCellMax) {</td>
      </tr>
      <tr>
        <td id="L1280" class="blob-num js-line-number" data-line-number="1280"></td>
        <td id="LC1280" class="blob-code js-file-line">                <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1281" class="blob-num js-line-number" data-line-number="1281"></td>
        <td id="LC1281" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1282" class="blob-num js-line-number" data-line-number="1282"></td>
        <td id="LC1282" class="blob-code js-file-line">            <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L1283" class="blob-num js-line-number" data-line-number="1283"></td>
        <td id="LC1283" class="blob-code js-file-line">                #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1284" class="blob-num js-line-number" data-line-number="1284"></td>
        <td id="LC1284" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>iGlobal = <span class="pl-c1">%d</span><span class="pl-cce">\n</span>iRang + 1 = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L1285" class="blob-num js-line-number" data-line-number="1285"></td>
        <td id="LC1285" class="blob-code js-file-line">                <span class="pl-s"><span class="pl-pds">&quot;</span>dbzRangeScale = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L1286" class="blob-num js-line-number" data-line-number="1286"></td>
        <td id="LC1286" class="blob-code js-file-line">                <span class="pl-s"><span class="pl-pds">&quot;</span>rCellMax = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L1287" class="blob-num js-line-number" data-line-number="1287"></td>
        <td id="LC1287" class="blob-code js-file-line">                <span class="pl-s"><span class="pl-pds">&quot;</span>(iRang + 1) * dbzRangeScale = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L1288" class="blob-num js-line-number" data-line-number="1288"></td>
        <td id="LC1288" class="blob-code js-file-line">                <span class="pl-s"><span class="pl-pds">&quot;</span>((iRang + 1) * dbzRangeScale &gt; rCellMax) = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L1289" class="blob-num js-line-number" data-line-number="1289"></td>
        <td id="LC1289" class="blob-code js-file-line">                <span class="pl-s"><span class="pl-pds">&quot;</span>dbg=<span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iGlobal,iRang + <span class="pl-c1">1</span>,dbzRangeScale,rCellMax,</td>
      </tr>
      <tr>
        <td id="L1290" class="blob-num js-line-number" data-line-number="1290"></td>
        <td id="LC1290" class="blob-code js-file-line">                (iRang + <span class="pl-c1">1</span>) * dbzRangeScale,</td>
      </tr>
      <tr>
        <td id="L1291" class="blob-num js-line-number" data-line-number="1291"></td>
        <td id="LC1291" class="blob-code js-file-line">                ((iRang + <span class="pl-c1">1</span>) * dbzRangeScale &gt; rCellMax),dbg);</td>
      </tr>
      <tr>
        <td id="L1292" class="blob-num js-line-number" data-line-number="1292"></td>
        <td id="LC1292" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L1293" class="blob-num js-line-number" data-line-number="1293"></td>
        <td id="LC1293" class="blob-code js-file-line">                dbg++;</td>
      </tr>
      <tr>
        <td id="L1294" class="blob-num js-line-number" data-line-number="1294"></td>
        <td id="LC1294" class="blob-code js-file-line">                #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1295" class="blob-num js-line-number" data-line-number="1295"></td>
        <td id="LC1295" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1296" class="blob-num js-line-number" data-line-number="1296"></td>
        <td id="LC1296" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1297" class="blob-num js-line-number" data-line-number="1297"></td>
        <td id="LC1297" class="blob-code js-file-line">            #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1298" class="blob-num js-line-number" data-line-number="1298"></td>
        <td id="LC1298" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>iGlobal = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iGlobal);</td>
      </tr>
      <tr>
        <td id="L1299" class="blob-num js-line-number" data-line-number="1299"></td>
        <td id="LC1299" class="blob-code js-file-line">            #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1300" class="blob-num js-line-number" data-line-number="1300"></td>
        <td id="LC1300" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1301" class="blob-num js-line-number" data-line-number="1301"></td>
        <td id="LC1301" class="blob-code js-file-line">            <span class="pl-k">if</span> (dbzImage[iGlobal] == (<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) dbzMissing) {</td>
      </tr>
      <tr>
        <td id="L1302" class="blob-num js-line-number" data-line-number="1302"></td>
        <td id="LC1302" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1303" class="blob-num js-line-number" data-line-number="1303"></td>
        <td id="LC1303" class="blob-code js-file-line">                #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1304" class="blob-num js-line-number" data-line-number="1304"></td>
        <td id="LC1304" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>dbzImage[<span class="pl-c1">%d</span>] == dbzMissing<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iGlobal);</td>
      </tr>
      <tr>
        <td id="L1305" class="blob-num js-line-number" data-line-number="1305"></td>
        <td id="LC1305" class="blob-code js-file-line">                #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1306" class="blob-num js-line-number" data-line-number="1306"></td>
        <td id="LC1306" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1307" class="blob-num js-line-number" data-line-number="1307"></td>
        <td id="LC1307" class="blob-code js-file-line">                <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1308" class="blob-num js-line-number" data-line-number="1308"></td>
        <td id="LC1308" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1309" class="blob-num js-line-number" data-line-number="1309"></td>
        <td id="LC1309" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1310" class="blob-num js-line-number" data-line-number="1310"></td>
        <td id="LC1310" class="blob-code js-file-line">            <span class="pl-k">if</span> (dbzImage[iGlobal] &lt; (<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) dbzThres) {</td>
      </tr>
      <tr>
        <td id="L1311" class="blob-num js-line-number" data-line-number="1311"></td>
        <td id="LC1311" class="blob-code js-file-line">                <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1312" class="blob-num js-line-number" data-line-number="1312"></td>
        <td id="LC1312" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1313" class="blob-num js-line-number" data-line-number="1313"></td>
        <td id="LC1313" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1314" class="blob-num js-line-number" data-line-number="1314"></td>
        <td id="LC1314" class="blob-code js-file-line">            <span class="pl-c">// count number of direct neighbors above threshold</span></td>
      </tr>
      <tr>
        <td id="L1315" class="blob-num js-line-number" data-line-number="1315"></td>
        <td id="LC1315" class="blob-code js-file-line">            count = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1316" class="blob-num js-line-number" data-line-number="1316"></td>
        <td id="LC1316" class="blob-code js-file-line">            <span class="pl-k">for</span> (iNeighborhood = <span class="pl-c1">0</span>; iNeighborhood &lt; nNeighborhood; iNeighborhood++) {</td>
      </tr>
      <tr>
        <td id="L1317" class="blob-num js-line-number" data-line-number="1317"></td>
        <td id="LC1317" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1318" class="blob-num js-line-number" data-line-number="1318"></td>
        <td id="LC1318" class="blob-code js-file-line">                iLocal = <span class="pl-c1">findNearbyGateIndex</span>(nAzim,nRang,iGlobal,nAzimNeighborhood,nRangNeighborhood,iNeighborhood);</td>
      </tr>
      <tr>
        <td id="L1319" class="blob-num js-line-number" data-line-number="1319"></td>
        <td id="LC1319" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1320" class="blob-num js-line-number" data-line-number="1320"></td>
        <td id="LC1320" class="blob-code js-file-line">                <span class="pl-k">if</span> (iLocal &lt; <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L1321" class="blob-num js-line-number" data-line-number="1321"></td>
        <td id="LC1321" class="blob-code js-file-line">                    <span class="pl-c">// iLocal below zero are error codes</span></td>
      </tr>
      <tr>
        <td id="L1322" class="blob-num js-line-number" data-line-number="1322"></td>
        <td id="LC1322" class="blob-code js-file-line">                    <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1323" class="blob-num js-line-number" data-line-number="1323"></td>
        <td id="LC1323" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1324" class="blob-num js-line-number" data-line-number="1324"></td>
        <td id="LC1324" class="blob-code js-file-line">                <span class="pl-k">if</span> (dbzImage[iLocal] &gt; dbzThres) {</td>
      </tr>
      <tr>
        <td id="L1325" class="blob-num js-line-number" data-line-number="1325"></td>
        <td id="LC1325" class="blob-code js-file-line">                    count++;</td>
      </tr>
      <tr>
        <td id="L1326" class="blob-num js-line-number" data-line-number="1326"></td>
        <td id="LC1326" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1327" class="blob-num js-line-number" data-line-number="1327"></td>
        <td id="LC1327" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1328" class="blob-num js-line-number" data-line-number="1328"></td>
        <td id="LC1328" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1329" class="blob-num js-line-number" data-line-number="1329"></td>
        <td id="LC1329" class="blob-code js-file-line">            <span class="pl-c">// when not enough qualified neighbors, continue</span></td>
      </tr>
      <tr>
        <td id="L1330" class="blob-num js-line-number" data-line-number="1330"></td>
        <td id="LC1330" class="blob-code js-file-line">            <span class="pl-k">if</span> (count - <span class="pl-c1">1</span> &lt; nNeighborsMin) {</td>
      </tr>
      <tr>
        <td id="L1331" class="blob-num js-line-number" data-line-number="1331"></td>
        <td id="LC1331" class="blob-code js-file-line">                <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1332" class="blob-num js-line-number" data-line-number="1332"></td>
        <td id="LC1332" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1333" class="blob-num js-line-number" data-line-number="1333"></td>
        <td id="LC1333" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1334" class="blob-num js-line-number" data-line-number="1334"></td>
        <td id="LC1334" class="blob-code js-file-line">            #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1335" class="blob-num js-line-number" data-line-number="1335"></td>
        <td id="LC1335" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>iGlobal = <span class="pl-c1">%d</span>, count = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iGlobal,count);</td>
      </tr>
      <tr>
        <td id="L1336" class="blob-num js-line-number" data-line-number="1336"></td>
        <td id="LC1336" class="blob-code js-file-line">            #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1337" class="blob-num js-line-number" data-line-number="1337"></td>
        <td id="LC1337" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1338" class="blob-num js-line-number" data-line-number="1338"></td>
        <td id="LC1338" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1339" class="blob-num js-line-number" data-line-number="1339"></td>
        <td id="LC1339" class="blob-code js-file-line">            <span class="pl-c">// Looking for horizontal, vertical, forward diagonal, and backward diagonal connections.</span></td>
      </tr>
      <tr>
        <td id="L1340" class="blob-num js-line-number" data-line-number="1340"></td>
        <td id="LC1340" class="blob-code js-file-line">            <span class="pl-k">for</span> (iNeighborhood = <span class="pl-c1">0</span>; iNeighborhood &lt; nHalfNeighborhood; iNeighborhood++) {</td>
      </tr>
      <tr>
        <td id="L1341" class="blob-num js-line-number" data-line-number="1341"></td>
        <td id="LC1341" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1342" class="blob-num js-line-number" data-line-number="1342"></td>
        <td id="LC1342" class="blob-code js-file-line">                iLocal = <span class="pl-c1">findNearbyGateIndex</span>(nAzim,nRang,iGlobal,nAzimNeighborhood,nRangNeighborhood,iNeighborhood);</td>
      </tr>
      <tr>
        <td id="L1343" class="blob-num js-line-number" data-line-number="1343"></td>
        <td id="LC1343" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1344" class="blob-num js-line-number" data-line-number="1344"></td>
        <td id="LC1344" class="blob-code js-file-line">                <span class="pl-k">if</span> (iLocal &lt; <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L1345" class="blob-num js-line-number" data-line-number="1345"></td>
        <td id="LC1345" class="blob-code js-file-line">                    <span class="pl-c">// iLocal less than zero are error codes</span></td>
      </tr>
      <tr>
        <td id="L1346" class="blob-num js-line-number" data-line-number="1346"></td>
        <td id="LC1346" class="blob-code js-file-line">                    <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1347" class="blob-num js-line-number" data-line-number="1347"></td>
        <td id="LC1347" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1348" class="blob-num js-line-number" data-line-number="1348"></td>
        <td id="LC1348" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1349" class="blob-num js-line-number" data-line-number="1349"></td>
        <td id="LC1349" class="blob-code js-file-line">                <span class="pl-c">// no connection found, go to next pixel within neighborhood</span></td>
      </tr>
      <tr>
        <td id="L1350" class="blob-num js-line-number" data-line-number="1350"></td>
        <td id="LC1350" class="blob-code js-file-line">                <span class="pl-k">if</span> (cellImage[iLocal] == cellImageInitialValue) {</td>
      </tr>
      <tr>
        <td id="L1351" class="blob-num js-line-number" data-line-number="1351"></td>
        <td id="LC1351" class="blob-code js-file-line">                    <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1352" class="blob-num js-line-number" data-line-number="1352"></td>
        <td id="LC1352" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1353" class="blob-num js-line-number" data-line-number="1353"></td>
        <td id="LC1353" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1354" class="blob-num js-line-number" data-line-number="1354"></td>
        <td id="LC1354" class="blob-code js-file-line">                <span class="pl-c">// if pixel still unassigned, assign same iCellIdentifier as connection</span></td>
      </tr>
      <tr>
        <td id="L1355" class="blob-num js-line-number" data-line-number="1355"></td>
        <td id="LC1355" class="blob-code js-file-line">                <span class="pl-k">if</span> (cellImage[iGlobal] == cellImageInitialValue) {</td>
      </tr>
      <tr>
        <td id="L1356" class="blob-num js-line-number" data-line-number="1356"></td>
        <td id="LC1356" class="blob-code js-file-line">                    cellImage[iGlobal] = cellImage[iLocal];</td>
      </tr>
      <tr>
        <td id="L1357" class="blob-num js-line-number" data-line-number="1357"></td>
        <td id="LC1357" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1358" class="blob-num js-line-number" data-line-number="1358"></td>
        <td id="LC1358" class="blob-code js-file-line">                <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L1359" class="blob-num js-line-number" data-line-number="1359"></td>
        <td id="LC1359" class="blob-code js-file-line">                    <span class="pl-c">// if connection found but pixel is already assigned a different iCellIdentifier:</span></td>
      </tr>
      <tr>
        <td id="L1360" class="blob-num js-line-number" data-line-number="1360"></td>
        <td id="LC1360" class="blob-code js-file-line">                    <span class="pl-k">if</span> (cellImage[iGlobal] != cellImage[iLocal]) {</td>
      </tr>
      <tr>
        <td id="L1361" class="blob-num js-line-number" data-line-number="1361"></td>
        <td id="LC1361" class="blob-code js-file-line">                        <span class="pl-c">// merging cells detected: replace all other occurences by value of connection:</span></td>
      </tr>
      <tr>
        <td id="L1362" class="blob-num js-line-number" data-line-number="1362"></td>
        <td id="LC1362" class="blob-code js-file-line">                        <span class="pl-k">for</span> (iGlobalOther = <span class="pl-c1">0</span>; iGlobalOther &lt; nGlobal; iGlobalOther++) {</td>
      </tr>
      <tr>
        <td id="L1363" class="blob-num js-line-number" data-line-number="1363"></td>
        <td id="LC1363" class="blob-code js-file-line">                            <span class="pl-k">if</span> (cellImage[iGlobalOther] == cellImage[iGlobal]) {</td>
      </tr>
      <tr>
        <td id="L1364" class="blob-num js-line-number" data-line-number="1364"></td>
        <td id="LC1364" class="blob-code js-file-line">                                cellImage[iGlobalOther] = cellImage[iLocal];</td>
      </tr>
      <tr>
        <td id="L1365" class="blob-num js-line-number" data-line-number="1365"></td>
        <td id="LC1365" class="blob-code js-file-line">                            }</td>
      </tr>
      <tr>
        <td id="L1366" class="blob-num js-line-number" data-line-number="1366"></td>
        <td id="LC1366" class="blob-code js-file-line">                            <span class="pl-c">// note: not all iCellIdentifier need to be used eventually</span></td>
      </tr>
      <tr>
        <td id="L1367" class="blob-num js-line-number" data-line-number="1367"></td>
        <td id="LC1367" class="blob-code js-file-line">                        }</td>
      </tr>
      <tr>
        <td id="L1368" class="blob-num js-line-number" data-line-number="1368"></td>
        <td id="LC1368" class="blob-code js-file-line">                    }</td>
      </tr>
      <tr>
        <td id="L1369" class="blob-num js-line-number" data-line-number="1369"></td>
        <td id="LC1369" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1370" class="blob-num js-line-number" data-line-number="1370"></td>
        <td id="LC1370" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1371" class="blob-num js-line-number" data-line-number="1371"></td>
        <td id="LC1371" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1372" class="blob-num js-line-number" data-line-number="1372"></td>
        <td id="LC1372" class="blob-code js-file-line">            <span class="pl-c">// When no connections are found, give a new number.</span></td>
      </tr>
      <tr>
        <td id="L1373" class="blob-num js-line-number" data-line-number="1373"></td>
        <td id="LC1373" class="blob-code js-file-line">            <span class="pl-k">if</span> (cellImage[iGlobal] == cellImageInitialValue) {</td>
      </tr>
      <tr>
        <td id="L1374" class="blob-num js-line-number" data-line-number="1374"></td>
        <td id="LC1374" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1375" class="blob-num js-line-number" data-line-number="1375"></td>
        <td id="LC1375" class="blob-code js-file-line">                #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1376" class="blob-num js-line-number" data-line-number="1376"></td>
        <td id="LC1376" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>new cell found...assigning number <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iCellIdentifier);</td>
      </tr>
      <tr>
        <td id="L1377" class="blob-num js-line-number" data-line-number="1377"></td>
        <td id="LC1377" class="blob-code js-file-line">                #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1378" class="blob-num js-line-number" data-line-number="1378"></td>
        <td id="LC1378" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1379" class="blob-num js-line-number" data-line-number="1379"></td>
        <td id="LC1379" class="blob-code js-file-line">                cellImage[iGlobal] = iCellIdentifier;</td>
      </tr>
      <tr>
        <td id="L1380" class="blob-num js-line-number" data-line-number="1380"></td>
        <td id="LC1380" class="blob-code js-file-line">                iCellIdentifier++;</td>
      </tr>
      <tr>
        <td id="L1381" class="blob-num js-line-number" data-line-number="1381"></td>
        <td id="LC1381" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1382" class="blob-num js-line-number" data-line-number="1382"></td>
        <td id="LC1382" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1383" class="blob-num js-line-number" data-line-number="1383"></td>
        <td id="LC1383" class="blob-code js-file-line">        } <span class="pl-c">// (iRang=0; iRang&lt;nRang; iRang++)</span></td>
      </tr>
      <tr>
        <td id="L1384" class="blob-num js-line-number" data-line-number="1384"></td>
        <td id="LC1384" class="blob-code js-file-line">    } <span class="pl-c">// for (iAzim=0; iAzim&lt;nAzim; iAzim++)</span></td>
      </tr>
      <tr>
        <td id="L1385" class="blob-num js-line-number" data-line-number="1385"></td>
        <td id="LC1385" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1386" class="blob-num js-line-number" data-line-number="1386"></td>
        <td id="LC1386" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1387" class="blob-num js-line-number" data-line-number="1387"></td>
        <td id="LC1387" class="blob-code js-file-line">    <span class="pl-c">// check whether a cell crosses the border of the array (remember that iAzim=0 is</span></td>
      </tr>
      <tr>
        <td id="L1388" class="blob-num js-line-number" data-line-number="1388"></td>
        <td id="LC1388" class="blob-code js-file-line">    <span class="pl-c">// adjacent to iAzim=nAzim-1):</span></td>
      </tr>
      <tr>
        <td id="L1389" class="blob-num js-line-number" data-line-number="1389"></td>
        <td id="LC1389" class="blob-code js-file-line">    iAzim = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1390" class="blob-num js-line-number" data-line-number="1390"></td>
        <td id="LC1390" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1391" class="blob-num js-line-number" data-line-number="1391"></td>
        <td id="LC1391" class="blob-code js-file-line">    <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L1392" class="blob-num js-line-number" data-line-number="1392"></td>
        <td id="LC1392" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1393" class="blob-num js-line-number" data-line-number="1393"></td>
        <td id="LC1393" class="blob-code js-file-line">        iGlobal = iRang + iAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L1394" class="blob-num js-line-number" data-line-number="1394"></td>
        <td id="LC1394" class="blob-code js-file-line">        <span class="pl-c">// index 1 in a 3x3 child array refers to the cell that is a direct neighbor of</span></td>
      </tr>
      <tr>
        <td id="L1395" class="blob-num js-line-number" data-line-number="1395"></td>
        <td id="LC1395" class="blob-code js-file-line">        <span class="pl-c">// iGlobal, but on the other side of the array (because the polar plot is wrapped</span></td>
      </tr>
      <tr>
        <td id="L1396" class="blob-num js-line-number" data-line-number="1396"></td>
        <td id="LC1396" class="blob-code js-file-line">        <span class="pl-c">// in the azimuth dimension):</span></td>
      </tr>
      <tr>
        <td id="L1397" class="blob-num js-line-number" data-line-number="1397"></td>
        <td id="LC1397" class="blob-code js-file-line">        iGlobalOther = <span class="pl-c1">findNearbyGateIndex</span>(nAzim,nRang,iGlobal,<span class="pl-c1">3</span>,<span class="pl-c1">3</span>,<span class="pl-c1">1</span>);</td>
      </tr>
      <tr>
        <td id="L1398" class="blob-num js-line-number" data-line-number="1398"></td>
        <td id="LC1398" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1399" class="blob-num js-line-number" data-line-number="1399"></td>
        <td id="LC1399" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1400" class="blob-num js-line-number" data-line-number="1400"></td>
        <td id="LC1400" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>iGlobal = <span class="pl-c1">%d</span>, iGlobalOther = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iGlobal,iGlobalOther);</td>
      </tr>
      <tr>
        <td id="L1401" class="blob-num js-line-number" data-line-number="1401"></td>
        <td id="LC1401" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1402" class="blob-num js-line-number" data-line-number="1402"></td>
        <td id="LC1402" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1403" class="blob-num js-line-number" data-line-number="1403"></td>
        <td id="LC1403" class="blob-code js-file-line">        cellIdentifierGlobal = cellImage[iGlobal];</td>
      </tr>
      <tr>
        <td id="L1404" class="blob-num js-line-number" data-line-number="1404"></td>
        <td id="LC1404" class="blob-code js-file-line">        cellIdentifierGlobalOther = cellImage[iGlobalOther];</td>
      </tr>
      <tr>
        <td id="L1405" class="blob-num js-line-number" data-line-number="1405"></td>
        <td id="LC1405" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellIdentifierGlobal != cellImageInitialValue &amp;&amp; cellIdentifierGlobalOther != cellImageInitialValue ) {</td>
      </tr>
      <tr>
        <td id="L1406" class="blob-num js-line-number" data-line-number="1406"></td>
        <td id="LC1406" class="blob-code js-file-line">            <span class="pl-c">// adjacent gates, both part of a cell -&gt; assign them the same identifier, i.e. assign</span></td>
      </tr>
      <tr>
        <td id="L1407" class="blob-num js-line-number" data-line-number="1407"></td>
        <td id="LC1407" class="blob-code js-file-line">            <span class="pl-c">// all elements of cellImage that are equal to cellImage[iGlobalOther] the value of</span></td>
      </tr>
      <tr>
        <td id="L1408" class="blob-num js-line-number" data-line-number="1408"></td>
        <td id="LC1408" class="blob-code js-file-line">            <span class="pl-c">// cellImage[iGlobal]</span></td>
      </tr>
      <tr>
        <td id="L1409" class="blob-num js-line-number" data-line-number="1409"></td>
        <td id="LC1409" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1410" class="blob-num js-line-number" data-line-number="1410"></td>
        <td id="LC1410" class="blob-code js-file-line">            <span class="pl-k">for</span> (iGlobalInner = <span class="pl-c1">0</span>; iGlobalInner &lt; nGlobal; iGlobalInner++) {</td>
      </tr>
      <tr>
        <td id="L1411" class="blob-num js-line-number" data-line-number="1411"></td>
        <td id="LC1411" class="blob-code js-file-line">                <span class="pl-k">if</span> (cellImage[iGlobalInner] == cellIdentifierGlobalOther) {</td>
      </tr>
      <tr>
        <td id="L1412" class="blob-num js-line-number" data-line-number="1412"></td>
        <td id="LC1412" class="blob-code js-file-line">                    cellImage[iGlobalInner] = cellIdentifierGlobal;</td>
      </tr>
      <tr>
        <td id="L1413" class="blob-num js-line-number" data-line-number="1413"></td>
        <td id="LC1413" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1414" class="blob-num js-line-number" data-line-number="1414"></td>
        <td id="LC1414" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1415" class="blob-num js-line-number" data-line-number="1415"></td>
        <td id="LC1415" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1416" class="blob-num js-line-number" data-line-number="1416"></td>
        <td id="LC1416" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1417" class="blob-num js-line-number" data-line-number="1417"></td>
        <td id="LC1417" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1418" class="blob-num js-line-number" data-line-number="1418"></td>
        <td id="LC1418" class="blob-code js-file-line">    <span class="pl-c">// Returning number of detected cells (including fringe/clutter)</span></td>
      </tr>
      <tr>
        <td id="L1419" class="blob-num js-line-number" data-line-number="1419"></td>
        <td id="LC1419" class="blob-code js-file-line">    nCells = iCellIdentifier;</td>
      </tr>
      <tr>
        <td id="L1420" class="blob-num js-line-number" data-line-number="1420"></td>
        <td id="LC1420" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1421" class="blob-num js-line-number" data-line-number="1421"></td>
        <td id="LC1421" class="blob-code js-file-line">    <span class="pl-k">return</span> nCells;</td>
      </tr>
      <tr>
        <td id="L1422" class="blob-num js-line-number" data-line-number="1422"></td>
        <td id="LC1422" class="blob-code js-file-line">} <span class="pl-c">// findWeatherCells</span></td>
      </tr>
      <tr>
        <td id="L1423" class="blob-num js-line-number" data-line-number="1423"></td>
        <td id="LC1423" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1424" class="blob-num js-line-number" data-line-number="1424"></td>
        <td id="LC1424" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1425" class="blob-num js-line-number" data-line-number="1425"></td>
        <td id="LC1425" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1426" class="blob-num js-line-number" data-line-number="1426"></td>
        <td id="LC1426" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1427" class="blob-num js-line-number" data-line-number="1427"></td>
        <td id="LC1427" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1428" class="blob-num js-line-number" data-line-number="1428"></td>
        <td id="LC1428" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1429" class="blob-num js-line-number" data-line-number="1429"></td>
        <td id="LC1429" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">findNearbyGateIndex</span>(<span class="pl-k">const</span> <span class="pl-k">int</span> nAzimParent, <span class="pl-k">const</span> <span class="pl-k">int</span> nRangParent, <span class="pl-k">const</span> <span class="pl-k">int</span> iParent,</td>
      </tr>
      <tr>
        <td id="L1430" class="blob-num js-line-number" data-line-number="1430"></td>
        <td id="LC1430" class="blob-code js-file-line">                        <span class="pl-k">const</span> <span class="pl-k">int</span> nAzimChild,  <span class="pl-k">const</span> <span class="pl-k">int</span> nRangChild,  <span class="pl-k">const</span> <span class="pl-k">int</span> iChild) {</td>
      </tr>
      <tr>
        <td id="L1431" class="blob-num js-line-number" data-line-number="1431"></td>
        <td id="LC1431" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1432" class="blob-num js-line-number" data-line-number="1432"></td>
        <td id="LC1432" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1433" class="blob-num js-line-number" data-line-number="1433"></td>
        <td id="LC1433" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1434" class="blob-num js-line-number" data-line-number="1434"></td>
        <td id="LC1434" class="blob-code js-file-line">    <span class="pl-k">if</span> (nRangChild%<span class="pl-c1">2</span> != <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L1435" class="blob-num js-line-number" data-line-number="1435"></td>
        <td id="LC1435" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1436" class="blob-num js-line-number" data-line-number="1436"></td>
        <td id="LC1436" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1437" class="blob-num js-line-number" data-line-number="1437"></td>
        <td id="LC1437" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>nRangChild must be an odd integer number.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1438" class="blob-num js-line-number" data-line-number="1438"></td>
        <td id="LC1438" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1439" class="blob-num js-line-number" data-line-number="1439"></td>
        <td id="LC1439" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1440" class="blob-num js-line-number" data-line-number="1440"></td>
        <td id="LC1440" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1441" class="blob-num js-line-number" data-line-number="1441"></td>
        <td id="LC1441" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1442" class="blob-num js-line-number" data-line-number="1442"></td>
        <td id="LC1442" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1443" class="blob-num js-line-number" data-line-number="1443"></td>
        <td id="LC1443" class="blob-code js-file-line">    <span class="pl-k">if</span> (nAzimChild%<span class="pl-c1">2</span> != <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L1444" class="blob-num js-line-number" data-line-number="1444"></td>
        <td id="LC1444" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1445" class="blob-num js-line-number" data-line-number="1445"></td>
        <td id="LC1445" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1446" class="blob-num js-line-number" data-line-number="1446"></td>
        <td id="LC1446" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>nAzimChild must be an odd integer number.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1447" class="blob-num js-line-number" data-line-number="1447"></td>
        <td id="LC1447" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1448" class="blob-num js-line-number" data-line-number="1448"></td>
        <td id="LC1448" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1449" class="blob-num js-line-number" data-line-number="1449"></td>
        <td id="LC1449" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">2</span>;</td>
      </tr>
      <tr>
        <td id="L1450" class="blob-num js-line-number" data-line-number="1450"></td>
        <td id="LC1450" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1451" class="blob-num js-line-number" data-line-number="1451"></td>
        <td id="LC1451" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1452" class="blob-num js-line-number" data-line-number="1452"></td>
        <td id="LC1452" class="blob-code js-file-line">    <span class="pl-k">if</span> (iChild &gt; nAzimChild * nRangChild - <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L1453" class="blob-num js-line-number" data-line-number="1453"></td>
        <td id="LC1453" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1454" class="blob-num js-line-number" data-line-number="1454"></td>
        <td id="LC1454" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1455" class="blob-num js-line-number" data-line-number="1455"></td>
        <td id="LC1455" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>iChild is outside the child window.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1456" class="blob-num js-line-number" data-line-number="1456"></td>
        <td id="LC1456" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1457" class="blob-num js-line-number" data-line-number="1457"></td>
        <td id="LC1457" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1458" class="blob-num js-line-number" data-line-number="1458"></td>
        <td id="LC1458" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">3</span>;</td>
      </tr>
      <tr>
        <td id="L1459" class="blob-num js-line-number" data-line-number="1459"></td>
        <td id="LC1459" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1460" class="blob-num js-line-number" data-line-number="1460"></td>
        <td id="LC1460" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1461" class="blob-num js-line-number" data-line-number="1461"></td>
        <td id="LC1461" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1462" class="blob-num js-line-number" data-line-number="1462"></td>
        <td id="LC1462" class="blob-code js-file-line">    <span class="pl-k">const</span> <span class="pl-k">int</span> iAzimParent = iParent / nRangParent;</td>
      </tr>
      <tr>
        <td id="L1463" class="blob-num js-line-number" data-line-number="1463"></td>
        <td id="LC1463" class="blob-code js-file-line">    <span class="pl-k">const</span> <span class="pl-k">int</span> iRangParent = iParent % nRangParent;</td>
      </tr>
      <tr>
        <td id="L1464" class="blob-num js-line-number" data-line-number="1464"></td>
        <td id="LC1464" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1465" class="blob-num js-line-number" data-line-number="1465"></td>
        <td id="LC1465" class="blob-code js-file-line">    <span class="pl-k">const</span> <span class="pl-k">int</span> iAzimChild = iChild / nRangChild;</td>
      </tr>
      <tr>
        <td id="L1466" class="blob-num js-line-number" data-line-number="1466"></td>
        <td id="LC1466" class="blob-code js-file-line">    <span class="pl-k">const</span> <span class="pl-k">int</span> iRangChild = iChild % nRangChild;</td>
      </tr>
      <tr>
        <td id="L1467" class="blob-num js-line-number" data-line-number="1467"></td>
        <td id="LC1467" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1468" class="blob-num js-line-number" data-line-number="1468"></td>
        <td id="LC1468" class="blob-code js-file-line">    <span class="pl-c">// the azimuth dimension is wrapped (polar plot); iAzim = 0 is adjacent to iAzim = nAzim-1:</span></td>
      </tr>
      <tr>
        <td id="L1469" class="blob-num js-line-number" data-line-number="1469"></td>
        <td id="LC1469" class="blob-code js-file-line">    <span class="pl-k">const</span> <span class="pl-k">int</span> iAzimReturn = (iAzimParent - nAzimChild/<span class="pl-c1">2</span> + iAzimChild + nAzimParent) % nAzimParent;</td>
      </tr>
      <tr>
        <td id="L1470" class="blob-num js-line-number" data-line-number="1470"></td>
        <td id="LC1470" class="blob-code js-file-line">    <span class="pl-k">const</span> <span class="pl-k">int</span> iRangReturn = iRangParent - nRangChild/<span class="pl-c1">2</span> + iRangChild;</td>
      </tr>
      <tr>
        <td id="L1471" class="blob-num js-line-number" data-line-number="1471"></td>
        <td id="LC1471" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1472" class="blob-num js-line-number" data-line-number="1472"></td>
        <td id="LC1472" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1473" class="blob-num js-line-number" data-line-number="1473"></td>
        <td id="LC1473" class="blob-code js-file-line">    <span class="pl-k">if</span> (iRangReturn &gt; nRangParent - <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L1474" class="blob-num js-line-number" data-line-number="1474"></td>
        <td id="LC1474" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1475" class="blob-num js-line-number" data-line-number="1475"></td>
        <td id="LC1475" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1476" class="blob-num js-line-number" data-line-number="1476"></td>
        <td id="LC1476" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>iChild is outside the parent array on the right-hand side.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1477" class="blob-num js-line-number" data-line-number="1477"></td>
        <td id="LC1477" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1478" class="blob-num js-line-number" data-line-number="1478"></td>
        <td id="LC1478" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1479" class="blob-num js-line-number" data-line-number="1479"></td>
        <td id="LC1479" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">4</span>;</td>
      </tr>
      <tr>
        <td id="L1480" class="blob-num js-line-number" data-line-number="1480"></td>
        <td id="LC1480" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1481" class="blob-num js-line-number" data-line-number="1481"></td>
        <td id="LC1481" class="blob-code js-file-line">    <span class="pl-k">if</span> (iRangReturn &lt; <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L1482" class="blob-num js-line-number" data-line-number="1482"></td>
        <td id="LC1482" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1483" class="blob-num js-line-number" data-line-number="1483"></td>
        <td id="LC1483" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1484" class="blob-num js-line-number" data-line-number="1484"></td>
        <td id="LC1484" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>iChild is outside the parent array on the left-hand side.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1485" class="blob-num js-line-number" data-line-number="1485"></td>
        <td id="LC1485" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1486" class="blob-num js-line-number" data-line-number="1486"></td>
        <td id="LC1486" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1487" class="blob-num js-line-number" data-line-number="1487"></td>
        <td id="LC1487" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">5</span>;</td>
      </tr>
      <tr>
        <td id="L1488" class="blob-num js-line-number" data-line-number="1488"></td>
        <td id="LC1488" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1489" class="blob-num js-line-number" data-line-number="1489"></td>
        <td id="LC1489" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1490" class="blob-num js-line-number" data-line-number="1490"></td>
        <td id="LC1490" class="blob-code js-file-line">    <span class="pl-k">return</span> iAzimReturn * nRangParent + iRangReturn;</td>
      </tr>
      <tr>
        <td id="L1491" class="blob-num js-line-number" data-line-number="1491"></td>
        <td id="LC1491" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1492" class="blob-num js-line-number" data-line-number="1492"></td>
        <td id="LC1492" class="blob-code js-file-line">} <span class="pl-c">// findNearbyGateIndex</span></td>
      </tr>
      <tr>
        <td id="L1493" class="blob-num js-line-number" data-line-number="1493"></td>
        <td id="LC1493" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1494" class="blob-num js-line-number" data-line-number="1494"></td>
        <td id="LC1494" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1495" class="blob-num js-line-number" data-line-number="1495"></td>
        <td id="LC1495" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1496" class="blob-num js-line-number" data-line-number="1496"></td>
        <td id="LC1496" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1497" class="blob-num js-line-number" data-line-number="1497"></td>
        <td id="LC1497" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1498" class="blob-num js-line-number" data-line-number="1498"></td>
        <td id="LC1498" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1499" class="blob-num js-line-number" data-line-number="1499"></td>
        <td id="LC1499" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1500" class="blob-num js-line-number" data-line-number="1500"></td>
        <td id="LC1500" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1501" class="blob-num js-line-number" data-line-number="1501"></td>
        <td id="LC1501" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">fringeCells</span>(<span class="pl-k">int</span> *cellImage, <span class="pl-k">int</span> nRang, <span class="pl-k">int</span> nAzim, <span class="pl-k">float</span> aScale, <span class="pl-k">float</span> rScale) {</td>
      </tr>
      <tr>
        <td id="L1502" class="blob-num js-line-number" data-line-number="1502"></td>
        <td id="LC1502" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1503" class="blob-num js-line-number" data-line-number="1503"></td>
        <td id="LC1503" class="blob-code js-file-line">    <span class="pl-c">// -------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1504" class="blob-num js-line-number" data-line-number="1504"></td>
        <td id="LC1504" class="blob-code js-file-line">    <span class="pl-c">// This function enlarges cells in cellImage by an additional fringe.         //</span></td>
      </tr>
      <tr>
        <td id="L1505" class="blob-num js-line-number" data-line-number="1505"></td>
        <td id="LC1505" class="blob-code js-file-line">    <span class="pl-c">// First a block around each pixel is searched for pixels within a distance   //</span></td>
      </tr>
      <tr>
        <td id="L1506" class="blob-num js-line-number" data-line-number="1506"></td>
        <td id="LC1506" class="blob-code js-file-line">    <span class="pl-c">// equal to &#39;fringeDist&#39;.                                                     //</span></td>
      </tr>
      <tr>
        <td id="L1507" class="blob-num js-line-number" data-line-number="1507"></td>
        <td id="LC1507" class="blob-code js-file-line">    <span class="pl-c">// -------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1508" class="blob-num js-line-number" data-line-number="1508"></td>
        <td id="LC1508" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1509" class="blob-num js-line-number" data-line-number="1509"></td>
        <td id="LC1509" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L1510" class="blob-num js-line-number" data-line-number="1510"></td>
        <td id="LC1510" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzim;</td>
      </tr>
      <tr>
        <td id="L1511" class="blob-num js-line-number" data-line-number="1511"></td>
        <td id="LC1511" class="blob-code js-file-line">    <span class="pl-k">int</span> iNeighborhood;</td>
      </tr>
      <tr>
        <td id="L1512" class="blob-num js-line-number" data-line-number="1512"></td>
        <td id="LC1512" class="blob-code js-file-line">    <span class="pl-k">int</span> nNeighborhood;</td>
      </tr>
      <tr>
        <td id="L1513" class="blob-num js-line-number" data-line-number="1513"></td>
        <td id="LC1513" class="blob-code js-file-line">    <span class="pl-k">int</span> iRangLocal;</td>
      </tr>
      <tr>
        <td id="L1514" class="blob-num js-line-number" data-line-number="1514"></td>
        <td id="LC1514" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzimLocal;</td>
      </tr>
      <tr>
        <td id="L1515" class="blob-num js-line-number" data-line-number="1515"></td>
        <td id="LC1515" class="blob-code js-file-line">    <span class="pl-k">int</span> iLocal;</td>
      </tr>
      <tr>
        <td id="L1516" class="blob-num js-line-number" data-line-number="1516"></td>
        <td id="LC1516" class="blob-code js-file-line">    <span class="pl-k">int</span> rBlock;</td>
      </tr>
      <tr>
        <td id="L1517" class="blob-num js-line-number" data-line-number="1517"></td>
        <td id="LC1517" class="blob-code js-file-line">    <span class="pl-k">int</span> aBlock;</td>
      </tr>
      <tr>
        <td id="L1518" class="blob-num js-line-number" data-line-number="1518"></td>
        <td id="LC1518" class="blob-code js-file-line">    <span class="pl-k">int</span> isEdge;</td>
      </tr>
      <tr>
        <td id="L1519" class="blob-num js-line-number" data-line-number="1519"></td>
        <td id="LC1519" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L1520" class="blob-num js-line-number" data-line-number="1520"></td>
        <td id="LC1520" class="blob-code js-file-line">    <span class="pl-k">float</span> theDist;</td>
      </tr>
      <tr>
        <td id="L1521" class="blob-num js-line-number" data-line-number="1521"></td>
        <td id="LC1521" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzimChild;</td>
      </tr>
      <tr>
        <td id="L1522" class="blob-num js-line-number" data-line-number="1522"></td>
        <td id="LC1522" class="blob-code js-file-line">    <span class="pl-k">int</span> nRangChild;</td>
      </tr>
      <tr>
        <td id="L1523" class="blob-num js-line-number" data-line-number="1523"></td>
        <td id="LC1523" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1524" class="blob-num js-line-number" data-line-number="1524"></td>
        <td id="LC1524" class="blob-code js-file-line">    <span class="pl-k">float</span> actualRange;</td>
      </tr>
      <tr>
        <td id="L1525" class="blob-num js-line-number" data-line-number="1525"></td>
        <td id="LC1525" class="blob-code js-file-line">    <span class="pl-k">float</span> circumferenceAtActualRange;</td>
      </tr>
      <tr>
        <td id="L1526" class="blob-num js-line-number" data-line-number="1526"></td>
        <td id="LC1526" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1527" class="blob-num js-line-number" data-line-number="1527"></td>
        <td id="LC1527" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1528" class="blob-num js-line-number" data-line-number="1528"></td>
        <td id="LC1528" class="blob-code js-file-line">    rBlock = <span class="pl-c1">ROUND</span>(fringeDist / rScale);</td>
      </tr>
      <tr>
        <td id="L1529" class="blob-num js-line-number" data-line-number="1529"></td>
        <td id="LC1529" class="blob-code js-file-line">    <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) {</td>
      </tr>
      <tr>
        <td id="L1530" class="blob-num js-line-number" data-line-number="1530"></td>
        <td id="LC1530" class="blob-code js-file-line">        <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L1531" class="blob-num js-line-number" data-line-number="1531"></td>
        <td id="LC1531" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1532" class="blob-num js-line-number" data-line-number="1532"></td>
        <td id="LC1532" class="blob-code js-file-line">            iGlobal = iRang + iAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L1533" class="blob-num js-line-number" data-line-number="1533"></td>
        <td id="LC1533" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1534" class="blob-num js-line-number" data-line-number="1534"></td>
        <td id="LC1534" class="blob-code js-file-line">            <span class="pl-k">if</span> (cellImage[iGlobal] &lt;= <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L1535" class="blob-num js-line-number" data-line-number="1535"></td>
        <td id="LC1535" class="blob-code js-file-line">                <span class="pl-k">continue</span>; <span class="pl-c">// with the next iGlobal; already fringe or not in cellImage</span></td>
      </tr>
      <tr>
        <td id="L1536" class="blob-num js-line-number" data-line-number="1536"></td>
        <td id="LC1536" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1537" class="blob-num js-line-number" data-line-number="1537"></td>
        <td id="LC1537" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1538" class="blob-num js-line-number" data-line-number="1538"></td>
        <td id="LC1538" class="blob-code js-file-line">            <span class="pl-c">// determine whether current pixel is a pixel on the edge of a cell</span></td>
      </tr>
      <tr>
        <td id="L1539" class="blob-num js-line-number" data-line-number="1539"></td>
        <td id="LC1539" class="blob-code js-file-line">            isEdge = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1540" class="blob-num js-line-number" data-line-number="1540"></td>
        <td id="LC1540" class="blob-code js-file-line">            <span class="pl-k">for</span> (iNeighborhood = <span class="pl-c1">0</span>; iNeighborhood &lt; <span class="pl-c1">9</span> &amp;&amp; isEdge == <span class="pl-c1">FALSE</span>; iNeighborhood++) {</td>
      </tr>
      <tr>
        <td id="L1541" class="blob-num js-line-number" data-line-number="1541"></td>
        <td id="LC1541" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1542" class="blob-num js-line-number" data-line-number="1542"></td>
        <td id="LC1542" class="blob-code js-file-line">                iLocal = <span class="pl-c1">findNearbyGateIndex</span>(nAzim,nRang,iGlobal,<span class="pl-c1">3</span>,<span class="pl-c1">3</span>,iNeighborhood);</td>
      </tr>
      <tr>
        <td id="L1543" class="blob-num js-line-number" data-line-number="1543"></td>
        <td id="LC1543" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1544" class="blob-num js-line-number" data-line-number="1544"></td>
        <td id="LC1544" class="blob-code js-file-line">                <span class="pl-k">if</span> (iLocal &lt; <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L1545" class="blob-num js-line-number" data-line-number="1545"></td>
        <td id="LC1545" class="blob-code js-file-line">                    <span class="pl-c">// iLocal less than zero are error codes</span></td>
      </tr>
      <tr>
        <td id="L1546" class="blob-num js-line-number" data-line-number="1546"></td>
        <td id="LC1546" class="blob-code js-file-line">                    <span class="pl-k">continue</span>; <span class="pl-c">// with the next iGlobal</span></td>
      </tr>
      <tr>
        <td id="L1547" class="blob-num js-line-number" data-line-number="1547"></td>
        <td id="LC1547" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1548" class="blob-num js-line-number" data-line-number="1548"></td>
        <td id="LC1548" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1549" class="blob-num js-line-number" data-line-number="1549"></td>
        <td id="LC1549" class="blob-code js-file-line">                <span class="pl-k">if</span> (cellImage[iLocal] &lt;= <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L1550" class="blob-num js-line-number" data-line-number="1550"></td>
        <td id="LC1550" class="blob-code js-file-line">                    isEdge = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1551" class="blob-num js-line-number" data-line-number="1551"></td>
        <td id="LC1551" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1552" class="blob-num js-line-number" data-line-number="1552"></td>
        <td id="LC1552" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1553" class="blob-num js-line-number" data-line-number="1553"></td>
        <td id="LC1553" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1554" class="blob-num js-line-number" data-line-number="1554"></td>
        <td id="LC1554" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L1555" class="blob-num js-line-number" data-line-number="1555"></td>
        <td id="LC1555" class="blob-code js-file-line">            <span class="pl-k">if</span> (isEdge == <span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L1556" class="blob-num js-line-number" data-line-number="1556"></td>
        <td id="LC1556" class="blob-code js-file-line">                <span class="pl-k">continue</span>; <span class="pl-c">// with the next iGlobal</span></td>
      </tr>
      <tr>
        <td id="L1557" class="blob-num js-line-number" data-line-number="1557"></td>
        <td id="LC1557" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L1558" class="blob-num js-line-number" data-line-number="1558"></td>
        <td id="LC1558" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1559" class="blob-num js-line-number" data-line-number="1559"></td>
        <td id="LC1559" class="blob-code js-file-line">            actualRange = (iRang+<span class="pl-c1">0.5</span>) * rScale;</td>
      </tr>
      <tr>
        <td id="L1560" class="blob-num js-line-number" data-line-number="1560"></td>
        <td id="LC1560" class="blob-code js-file-line">            circumferenceAtActualRange = <span class="pl-c1">2</span> * PI * actualRange;</td>
      </tr>
      <tr>
        <td id="L1561" class="blob-num js-line-number" data-line-number="1561"></td>
        <td id="LC1561" class="blob-code js-file-line">            aBlock = (fringeDist / circumferenceAtActualRange) * nAzim;</td>
      </tr>
      <tr>
        <td id="L1562" class="blob-num js-line-number" data-line-number="1562"></td>
        <td id="LC1562" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1563" class="blob-num js-line-number" data-line-number="1563"></td>
        <td id="LC1563" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1564" class="blob-num js-line-number" data-line-number="1564"></td>
        <td id="LC1564" class="blob-code js-file-line">            #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L1565" class="blob-num js-line-number" data-line-number="1565"></td>
        <td id="LC1565" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>actualRange = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,actualRange);</td>
      </tr>
      <tr>
        <td id="L1566" class="blob-num js-line-number" data-line-number="1566"></td>
        <td id="LC1566" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>circumferenceAtActualRange = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,circumferenceAtActualRange);</td>
      </tr>
      <tr>
        <td id="L1567" class="blob-num js-line-number" data-line-number="1567"></td>
        <td id="LC1567" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>fringeDist / circumferenceAtActualRange = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,fringeDist / circumferenceAtActualRange);</td>
      </tr>
      <tr>
        <td id="L1568" class="blob-num js-line-number" data-line-number="1568"></td>
        <td id="LC1568" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>aBlock = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>, aBlock);</td>
      </tr>
      <tr>
        <td id="L1569" class="blob-num js-line-number" data-line-number="1569"></td>
        <td id="LC1569" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>rBlock = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>, rBlock);</td>
      </tr>
      <tr>
        <td id="L1570" class="blob-num js-line-number" data-line-number="1570"></td>
        <td id="LC1570" class="blob-code js-file-line">            #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L1571" class="blob-num js-line-number" data-line-number="1571"></td>
        <td id="LC1571" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1572" class="blob-num js-line-number" data-line-number="1572"></td>
        <td id="LC1572" class="blob-code js-file-line">            nAzimChild = <span class="pl-c1">2</span> * aBlock + <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1573" class="blob-num js-line-number" data-line-number="1573"></td>
        <td id="LC1573" class="blob-code js-file-line">            nRangChild = <span class="pl-c1">2</span> * rBlock + <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1574" class="blob-num js-line-number" data-line-number="1574"></td>
        <td id="LC1574" class="blob-code js-file-line">            nNeighborhood = nAzimChild * nRangChild;</td>
      </tr>
      <tr>
        <td id="L1575" class="blob-num js-line-number" data-line-number="1575"></td>
        <td id="LC1575" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1576" class="blob-num js-line-number" data-line-number="1576"></td>
        <td id="LC1576" class="blob-code js-file-line">            <span class="pl-k">for</span> (iNeighborhood = <span class="pl-c1">0</span>; iNeighborhood &lt; nNeighborhood; iNeighborhood++) {</td>
      </tr>
      <tr>
        <td id="L1577" class="blob-num js-line-number" data-line-number="1577"></td>
        <td id="LC1577" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1578" class="blob-num js-line-number" data-line-number="1578"></td>
        <td id="LC1578" class="blob-code js-file-line">                iLocal = <span class="pl-c1">findNearbyGateIndex</span>(nAzim,nRang,iGlobal,nAzimChild,nRangChild,iNeighborhood);</td>
      </tr>
      <tr>
        <td id="L1579" class="blob-num js-line-number" data-line-number="1579"></td>
        <td id="LC1579" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1580" class="blob-num js-line-number" data-line-number="1580"></td>
        <td id="LC1580" class="blob-code js-file-line">                <span class="pl-k">if</span> (iLocal &lt; <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L1581" class="blob-num js-line-number" data-line-number="1581"></td>
        <td id="LC1581" class="blob-code js-file-line">                    <span class="pl-c">// iLocal less than zero are error codes</span></td>
      </tr>
      <tr>
        <td id="L1582" class="blob-num js-line-number" data-line-number="1582"></td>
        <td id="LC1582" class="blob-code js-file-line">                    <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1583" class="blob-num js-line-number" data-line-number="1583"></td>
        <td id="LC1583" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1584" class="blob-num js-line-number" data-line-number="1584"></td>
        <td id="LC1584" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1585" class="blob-num js-line-number" data-line-number="1585"></td>
        <td id="LC1585" class="blob-code js-file-line">                iAzimLocal = iLocal / nRang;</td>
      </tr>
      <tr>
        <td id="L1586" class="blob-num js-line-number" data-line-number="1586"></td>
        <td id="LC1586" class="blob-code js-file-line">                iRangLocal = iLocal % nRang;</td>
      </tr>
      <tr>
        <td id="L1587" class="blob-num js-line-number" data-line-number="1587"></td>
        <td id="LC1587" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1588" class="blob-num js-line-number" data-line-number="1588"></td>
        <td id="LC1588" class="blob-code js-file-line">                <span class="pl-c">// if not within range or already in cellImage or already a fringe, do nothing</span></td>
      </tr>
      <tr>
        <td id="L1589" class="blob-num js-line-number" data-line-number="1589"></td>
        <td id="LC1589" class="blob-code js-file-line">                theDist = <span class="pl-c1">calcDist</span>(iRang, iAzim, iRangLocal, iAzimLocal, rScale, aScale);</td>
      </tr>
      <tr>
        <td id="L1590" class="blob-num js-line-number" data-line-number="1590"></td>
        <td id="LC1590" class="blob-code js-file-line">                <span class="pl-k">if</span> (theDist &gt; fringeDist || cellImage[iLocal] &gt;= <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L1591" class="blob-num js-line-number" data-line-number="1591"></td>
        <td id="LC1591" class="blob-code js-file-line">                    <span class="pl-k">continue</span>; <span class="pl-c">// with the next iGlobal</span></td>
      </tr>
      <tr>
        <td id="L1592" class="blob-num js-line-number" data-line-number="1592"></td>
        <td id="LC1592" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L1593" class="blob-num js-line-number" data-line-number="1593"></td>
        <td id="LC1593" class="blob-code js-file-line">                <span class="pl-c">// include pixel (iRangLocal,iAzimLocal) in fringe</span></td>
      </tr>
      <tr>
        <td id="L1594" class="blob-num js-line-number" data-line-number="1594"></td>
        <td id="LC1594" class="blob-code js-file-line">                cellImage[iLocal] = <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1595" class="blob-num js-line-number" data-line-number="1595"></td>
        <td id="LC1595" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L1596" class="blob-num js-line-number" data-line-number="1596"></td>
        <td id="LC1596" class="blob-code js-file-line">            } <span class="pl-c">// (iNeighborhood = 0; iNeighborhood &lt; nNeighborhood; iNeighborhood++)</span></td>
      </tr>
      <tr>
        <td id="L1597" class="blob-num js-line-number" data-line-number="1597"></td>
        <td id="LC1597" class="blob-code js-file-line">        } <span class="pl-c">// (iRang = 0; iRang &lt; nRang; iRang++)</span></td>
      </tr>
      <tr>
        <td id="L1598" class="blob-num js-line-number" data-line-number="1598"></td>
        <td id="LC1598" class="blob-code js-file-line">    } <span class="pl-c">// (iAzim = 0; iAzim &lt; nAzim; iAzim++)</span></td>
      </tr>
      <tr>
        <td id="L1599" class="blob-num js-line-number" data-line-number="1599"></td>
        <td id="LC1599" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1600" class="blob-num js-line-number" data-line-number="1600"></td>
        <td id="LC1600" class="blob-code js-file-line">    <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L1601" class="blob-num js-line-number" data-line-number="1601"></td>
        <td id="LC1601" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1602" class="blob-num js-line-number" data-line-number="1602"></td>
        <td id="LC1602" class="blob-code js-file-line">} <span class="pl-c">// fringeCells</span></td>
      </tr>
      <tr>
        <td id="L1603" class="blob-num js-line-number" data-line-number="1603"></td>
        <td id="LC1603" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1604" class="blob-num js-line-number" data-line-number="1604"></td>
        <td id="LC1604" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1605" class="blob-num js-line-number" data-line-number="1605"></td>
        <td id="LC1605" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1606" class="blob-num js-line-number" data-line-number="1606"></td>
        <td id="LC1606" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1607" class="blob-num js-line-number" data-line-number="1607"></td>
        <td id="LC1607" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">getListOfSelectedGates</span>(<span class="pl-k">const</span> SCANMETA* vradMeta, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *vradImage, </td>
      </tr>
      <tr>
        <td id="L1608" class="blob-num js-line-number" data-line-number="1608"></td>
        <td id="LC1608" class="blob-code js-file-line">                           <span class="pl-k">const</span> SCANMETA* dbzMeta, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span> *dbzImage,</td>
      </tr>
      <tr>
        <td id="L1609" class="blob-num js-line-number" data-line-number="1609"></td>
        <td id="LC1609" class="blob-code js-file-line">                           <span class="pl-k">const</span> <span class="pl-k">int</span> *cellImage,</td>
      </tr>
      <tr>
        <td id="L1610" class="blob-num js-line-number" data-line-number="1610"></td>
        <td id="LC1610" class="blob-code js-file-line">                           <span class="pl-k">const</span> <span class="pl-k">float</span> altitudeMin, <span class="pl-k">const</span> <span class="pl-k">float</span> altitudeMax,</td>
      </tr>
      <tr>
        <td id="L1611" class="blob-num js-line-number" data-line-number="1611"></td>
        <td id="LC1611" class="blob-code js-file-line">                           <span class="pl-k">float</span>* points, <span class="pl-k">int</span> iRowPoints) {</td>
      </tr>
      <tr>
        <td id="L1612" class="blob-num js-line-number" data-line-number="1612"></td>
        <td id="LC1612" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1613" class="blob-num js-line-number" data-line-number="1613"></td>
        <td id="LC1613" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1614" class="blob-num js-line-number" data-line-number="1614"></td>
        <td id="LC1614" class="blob-code js-file-line">    <span class="pl-c">// Write combinations of an azimuth angle, an elevation angle, an      // </span></td>
      </tr>
      <tr>
        <td id="L1615" class="blob-num js-line-number" data-line-number="1615"></td>
        <td id="LC1615" class="blob-code js-file-line">    <span class="pl-c">// observed vrad value, an observed dbz value, and a cell identifier   //</span></td>
      </tr>
      <tr>
        <td id="L1616" class="blob-num js-line-number" data-line-number="1616"></td>
        <td id="LC1616" class="blob-code js-file-line">    <span class="pl-c">// value into an external larger list.                                 //</span></td>
      </tr>
      <tr>
        <td id="L1617" class="blob-num js-line-number" data-line-number="1617"></td>
        <td id="LC1617" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L1618" class="blob-num js-line-number" data-line-number="1618"></td>
        <td id="LC1618" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1619" class="blob-num js-line-number" data-line-number="1619"></td>
        <td id="LC1619" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzim;</td>
      </tr>
      <tr>
        <td id="L1620" class="blob-num js-line-number" data-line-number="1620"></td>
        <td id="LC1620" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L1621" class="blob-num js-line-number" data-line-number="1621"></td>
        <td id="LC1621" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L1622" class="blob-num js-line-number" data-line-number="1622"></td>
        <td id="LC1622" class="blob-code js-file-line">    <span class="pl-k">int</span> nRang;</td>
      </tr>
      <tr>
        <td id="L1623" class="blob-num js-line-number" data-line-number="1623"></td>
        <td id="LC1623" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzim;</td>
      </tr>
      <tr>
        <td id="L1624" class="blob-num js-line-number" data-line-number="1624"></td>
        <td id="LC1624" class="blob-code js-file-line">    <span class="pl-k">int</span> nPointsWritten;</td>
      </tr>
      <tr>
        <td id="L1625" class="blob-num js-line-number" data-line-number="1625"></td>
        <td id="LC1625" class="blob-code js-file-line">    <span class="pl-k">const</span> <span class="pl-k">int</span> nColsPoints = <span class="pl-c1">6</span>;</td>
      </tr>
      <tr>
        <td id="L1626" class="blob-num js-line-number" data-line-number="1626"></td>
        <td id="LC1626" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1627" class="blob-num js-line-number" data-line-number="1627"></td>
        <td id="LC1627" class="blob-code js-file-line">    <span class="pl-k">float</span> gateHeight;</td>
      </tr>
      <tr>
        <td id="L1628" class="blob-num js-line-number" data-line-number="1628"></td>
        <td id="LC1628" class="blob-code js-file-line">    <span class="pl-k">float</span> gateRange;</td>
      </tr>
      <tr>
        <td id="L1629" class="blob-num js-line-number" data-line-number="1629"></td>
        <td id="LC1629" class="blob-code js-file-line">    <span class="pl-k">float</span> gateAzim;</td>
      </tr>
      <tr>
        <td id="L1630" class="blob-num js-line-number" data-line-number="1630"></td>
        <td id="LC1630" class="blob-code js-file-line">    <span class="pl-k">float</span> rangeScale;</td>
      </tr>
      <tr>
        <td id="L1631" class="blob-num js-line-number" data-line-number="1631"></td>
        <td id="LC1631" class="blob-code js-file-line">    <span class="pl-k">float</span> azimuthScale;</td>
      </tr>
      <tr>
        <td id="L1632" class="blob-num js-line-number" data-line-number="1632"></td>
        <td id="LC1632" class="blob-code js-file-line">    <span class="pl-k">float</span> elevAngle;</td>
      </tr>
      <tr>
        <td id="L1633" class="blob-num js-line-number" data-line-number="1633"></td>
        <td id="LC1633" class="blob-code js-file-line">    <span class="pl-k">float</span> radarHeight;</td>
      </tr>
      <tr>
        <td id="L1634" class="blob-num js-line-number" data-line-number="1634"></td>
        <td id="LC1634" class="blob-code js-file-line">    <span class="pl-k">float</span> vradValueOffset;</td>
      </tr>
      <tr>
        <td id="L1635" class="blob-num js-line-number" data-line-number="1635"></td>
        <td id="LC1635" class="blob-code js-file-line">    <span class="pl-k">float</span> vradValueScale;</td>
      </tr>
      <tr>
        <td id="L1636" class="blob-num js-line-number" data-line-number="1636"></td>
        <td id="LC1636" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzValueOffset;</td>
      </tr>
      <tr>
        <td id="L1637" class="blob-num js-line-number" data-line-number="1637"></td>
        <td id="LC1637" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzValueScale;</td>
      </tr>
      <tr>
        <td id="L1638" class="blob-num js-line-number" data-line-number="1638"></td>
        <td id="LC1638" class="blob-code js-file-line">    <span class="pl-k">float</span> vradValue;</td>
      </tr>
      <tr>
        <td id="L1639" class="blob-num js-line-number" data-line-number="1639"></td>
        <td id="LC1639" class="blob-code js-file-line">    <span class="pl-k">float</span> dbzValue;</td>
      </tr>
      <tr>
        <td id="L1640" class="blob-num js-line-number" data-line-number="1640"></td>
        <td id="LC1640" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1641" class="blob-num js-line-number" data-line-number="1641"></td>
        <td id="LC1641" class="blob-code js-file-line">    nPointsWritten = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1642" class="blob-num js-line-number" data-line-number="1642"></td>
        <td id="LC1642" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1643" class="blob-num js-line-number" data-line-number="1643"></td>
        <td id="LC1643" class="blob-code js-file-line">    nRang = vradMeta-&gt;nRang;</td>
      </tr>
      <tr>
        <td id="L1644" class="blob-num js-line-number" data-line-number="1644"></td>
        <td id="LC1644" class="blob-code js-file-line">    nAzim = vradMeta-&gt;nAzim;</td>
      </tr>
      <tr>
        <td id="L1645" class="blob-num js-line-number" data-line-number="1645"></td>
        <td id="LC1645" class="blob-code js-file-line">    rangeScale = vradMeta-&gt;rangeScale;</td>
      </tr>
      <tr>
        <td id="L1646" class="blob-num js-line-number" data-line-number="1646"></td>
        <td id="LC1646" class="blob-code js-file-line">    azimuthScale = vradMeta-&gt;azimScale;</td>
      </tr>
      <tr>
        <td id="L1647" class="blob-num js-line-number" data-line-number="1647"></td>
        <td id="LC1647" class="blob-code js-file-line">    elevAngle = vradMeta-&gt;elev;</td>
      </tr>
      <tr>
        <td id="L1648" class="blob-num js-line-number" data-line-number="1648"></td>
        <td id="LC1648" class="blob-code js-file-line">    radarHeight = vradMeta-&gt;heig;</td>
      </tr>
      <tr>
        <td id="L1649" class="blob-num js-line-number" data-line-number="1649"></td>
        <td id="LC1649" class="blob-code js-file-line">    vradValueOffset = vradMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L1650" class="blob-num js-line-number" data-line-number="1650"></td>
        <td id="LC1650" class="blob-code js-file-line">    vradValueScale = vradMeta-&gt;valueScale;</td>
      </tr>
      <tr>
        <td id="L1651" class="blob-num js-line-number" data-line-number="1651"></td>
        <td id="LC1651" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1652" class="blob-num js-line-number" data-line-number="1652"></td>
        <td id="LC1652" class="blob-code js-file-line">    dbzValueOffset = dbzMeta-&gt;valueOffset;</td>
      </tr>
      <tr>
        <td id="L1653" class="blob-num js-line-number" data-line-number="1653"></td>
        <td id="LC1653" class="blob-code js-file-line">    dbzValueScale = dbzMeta-&gt;valueScale;</td>
      </tr>
      <tr>
        <td id="L1654" class="blob-num js-line-number" data-line-number="1654"></td>
        <td id="LC1654" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1655" class="blob-num js-line-number" data-line-number="1655"></td>
        <td id="LC1655" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1656" class="blob-num js-line-number" data-line-number="1656"></td>
        <td id="LC1656" class="blob-code js-file-line">    <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L1657" class="blob-num js-line-number" data-line-number="1657"></td>
        <td id="LC1657" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1658" class="blob-num js-line-number" data-line-number="1658"></td>
        <td id="LC1658" class="blob-code js-file-line">        <span class="pl-c">// so gateRange represents a distance along the view direction (not necessarily horizontal)</span></td>
      </tr>
      <tr>
        <td id="L1659" class="blob-num js-line-number" data-line-number="1659"></td>
        <td id="LC1659" class="blob-code js-file-line">        gateRange = ((<span class="pl-k">float</span>) iRang + <span class="pl-c1">0</span>.<span class="pl-c1">5f</span>) * rangeScale;</td>
      </tr>
      <tr>
        <td id="L1660" class="blob-num js-line-number" data-line-number="1660"></td>
        <td id="LC1660" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1661" class="blob-num js-line-number" data-line-number="1661"></td>
        <td id="LC1661" class="blob-code js-file-line">        <span class="pl-c">// note that &quot;sin(elevAngle*DEG2RAD)&quot; is equivalent to = &quot;cos((90 - elevAngle)*DEG2RAD)&quot;:</span></td>
      </tr>
      <tr>
        <td id="L1662" class="blob-num js-line-number" data-line-number="1662"></td>
        <td id="LC1662" class="blob-code js-file-line">        gateHeight = gateRange * (<span class="pl-k">float</span>) <span class="pl-c1">sin</span>(elevAngle*DEG2RAD) + radarHeight;</td>
      </tr>
      <tr>
        <td id="L1663" class="blob-num js-line-number" data-line-number="1663"></td>
        <td id="LC1663" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1664" class="blob-num js-line-number" data-line-number="1664"></td>
        <td id="LC1664" class="blob-code js-file-line">        <span class="pl-k">if</span> (gateRange &lt; rangeMin || gateRange &gt; rangeMax) {</td>
      </tr>
      <tr>
        <td id="L1665" class="blob-num js-line-number" data-line-number="1665"></td>
        <td id="LC1665" class="blob-code js-file-line">            <span class="pl-c">// the current gate is either </span></td>
      </tr>
      <tr>
        <td id="L1666" class="blob-num js-line-number" data-line-number="1666"></td>
        <td id="LC1666" class="blob-code js-file-line">            <span class="pl-c">// (1) too close to the radar; or</span></td>
      </tr>
      <tr>
        <td id="L1667" class="blob-num js-line-number" data-line-number="1667"></td>
        <td id="LC1667" class="blob-code js-file-line">            <span class="pl-c">// (2) too far away.</span></td>
      </tr>
      <tr>
        <td id="L1668" class="blob-num js-line-number" data-line-number="1668"></td>
        <td id="LC1668" class="blob-code js-file-line">            <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1669" class="blob-num js-line-number" data-line-number="1669"></td>
        <td id="LC1669" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1670" class="blob-num js-line-number" data-line-number="1670"></td>
        <td id="LC1670" class="blob-code js-file-line">        <span class="pl-k">if</span> (gateHeight &lt; altitudeMin || gateHeight &gt; altitudeMax) {</td>
      </tr>
      <tr>
        <td id="L1671" class="blob-num js-line-number" data-line-number="1671"></td>
        <td id="LC1671" class="blob-code js-file-line">            <span class="pl-c">// if the height of the middle of the current gate is too far away from</span></td>
      </tr>
      <tr>
        <td id="L1672" class="blob-num js-line-number" data-line-number="1672"></td>
        <td id="LC1672" class="blob-code js-file-line">            <span class="pl-c">// the requested height, continue with the next gate</span></td>
      </tr>
      <tr>
        <td id="L1673" class="blob-num js-line-number" data-line-number="1673"></td>
        <td id="LC1673" class="blob-code js-file-line">            <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L1674" class="blob-num js-line-number" data-line-number="1674"></td>
        <td id="LC1674" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1675" class="blob-num js-line-number" data-line-number="1675"></td>
        <td id="LC1675" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1676" class="blob-num js-line-number" data-line-number="1676"></td>
        <td id="LC1676" class="blob-code js-file-line">        <span class="pl-c">// the gates at this range and elevation angle are within bounds,</span></td>
      </tr>
      <tr>
        <td id="L1677" class="blob-num js-line-number" data-line-number="1677"></td>
        <td id="LC1677" class="blob-code js-file-line">        <span class="pl-c">// include their data in the &#39;points&#39; array:</span></td>
      </tr>
      <tr>
        <td id="L1678" class="blob-num js-line-number" data-line-number="1678"></td>
        <td id="LC1678" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1679" class="blob-num js-line-number" data-line-number="1679"></td>
        <td id="LC1679" class="blob-code js-file-line">        <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) {</td>
      </tr>
      <tr>
        <td id="L1680" class="blob-num js-line-number" data-line-number="1680"></td>
        <td id="LC1680" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1681" class="blob-num js-line-number" data-line-number="1681"></td>
        <td id="LC1681" class="blob-code js-file-line">            iGlobal = iRang + iAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L1682" class="blob-num js-line-number" data-line-number="1682"></td>
        <td id="LC1682" class="blob-code js-file-line">            gateAzim = ((<span class="pl-k">float</span>) iAzim + <span class="pl-c1">0</span>.<span class="pl-c1">5f</span>) * azimuthScale;</td>
      </tr>
      <tr>
        <td id="L1683" class="blob-num js-line-number" data-line-number="1683"></td>
        <td id="LC1683" class="blob-code js-file-line">            vradValue = vradValueScale * (<span class="pl-k">float</span>) vradImage[iGlobal] + vradValueOffset;</td>
      </tr>
      <tr>
        <td id="L1684" class="blob-num js-line-number" data-line-number="1684"></td>
        <td id="LC1684" class="blob-code js-file-line">            dbzValue = dbzValueScale * (<span class="pl-k">float</span>) dbzImage[iGlobal] + dbzValueOffset;</td>
      </tr>
      <tr>
        <td id="L1685" class="blob-num js-line-number" data-line-number="1685"></td>
        <td id="LC1685" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1686" class="blob-num js-line-number" data-line-number="1686"></td>
        <td id="LC1686" class="blob-code js-file-line">            <span class="pl-c">// store the location as an azimuth angle, elevation angle combination</span></td>
      </tr>
      <tr>
        <td id="L1687" class="blob-num js-line-number" data-line-number="1687"></td>
        <td id="LC1687" class="blob-code js-file-line">            points[iRowPoints * nColsPoints + <span class="pl-c1">0</span>] = gateAzim;</td>
      </tr>
      <tr>
        <td id="L1688" class="blob-num js-line-number" data-line-number="1688"></td>
        <td id="LC1688" class="blob-code js-file-line">            points[iRowPoints * nColsPoints + <span class="pl-c1">1</span>] = elevAngle;</td>
      </tr>
      <tr>
        <td id="L1689" class="blob-num js-line-number" data-line-number="1689"></td>
        <td id="LC1689" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1690" class="blob-num js-line-number" data-line-number="1690"></td>
        <td id="LC1690" class="blob-code js-file-line">            <span class="pl-c">// also store the dbz value --useful when estimating the bird density</span></td>
      </tr>
      <tr>
        <td id="L1691" class="blob-num js-line-number" data-line-number="1691"></td>
        <td id="LC1691" class="blob-code js-file-line">            points[iRowPoints * nColsPoints + <span class="pl-c1">2</span>] = dbzValue;</td>
      </tr>
      <tr>
        <td id="L1692" class="blob-num js-line-number" data-line-number="1692"></td>
        <td id="LC1692" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L1693" class="blob-num js-line-number" data-line-number="1693"></td>
        <td id="LC1693" class="blob-code js-file-line">            <span class="pl-c">// store the corresponding observed vrad value</span></td>
      </tr>
      <tr>
        <td id="L1694" class="blob-num js-line-number" data-line-number="1694"></td>
        <td id="LC1694" class="blob-code js-file-line">            points[iRowPoints * nColsPoints + <span class="pl-c1">3</span>] = vradValue;</td>
      </tr>
      <tr>
        <td id="L1695" class="blob-num js-line-number" data-line-number="1695"></td>
        <td id="LC1695" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1696" class="blob-num js-line-number" data-line-number="1696"></td>
        <td id="LC1696" class="blob-code js-file-line">            <span class="pl-c">// store the corresponding cellImage value</span></td>
      </tr>
      <tr>
        <td id="L1697" class="blob-num js-line-number" data-line-number="1697"></td>
        <td id="LC1697" class="blob-code js-file-line">            points[iRowPoints * nColsPoints + <span class="pl-c1">4</span>] = (<span class="pl-k">float</span>) cellImage[iGlobal];</td>
      </tr>
      <tr>
        <td id="L1698" class="blob-num js-line-number" data-line-number="1698"></td>
        <td id="LC1698" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1699" class="blob-num js-line-number" data-line-number="1699"></td>
        <td id="LC1699" class="blob-code js-file-line">            <span class="pl-c">// set the gateCode to zero for now</span></td>
      </tr>
      <tr>
        <td id="L1700" class="blob-num js-line-number" data-line-number="1700"></td>
        <td id="LC1700" class="blob-code js-file-line">            points[iRowPoints * nColsPoints + <span class="pl-c1">5</span>] = (<span class="pl-k">float</span>) <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1701" class="blob-num js-line-number" data-line-number="1701"></td>
        <td id="LC1701" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1702" class="blob-num js-line-number" data-line-number="1702"></td>
        <td id="LC1702" class="blob-code js-file-line">            <span class="pl-c">// raise the row counter by 1</span></td>
      </tr>
      <tr>
        <td id="L1703" class="blob-num js-line-number" data-line-number="1703"></td>
        <td id="LC1703" class="blob-code js-file-line">            iRowPoints += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1704" class="blob-num js-line-number" data-line-number="1704"></td>
        <td id="LC1704" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L1705" class="blob-num js-line-number" data-line-number="1705"></td>
        <td id="LC1705" class="blob-code js-file-line">            <span class="pl-c">// raise number of points written by 1</span></td>
      </tr>
      <tr>
        <td id="L1706" class="blob-num js-line-number" data-line-number="1706"></td>
        <td id="LC1706" class="blob-code js-file-line">            nPointsWritten += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1707" class="blob-num js-line-number" data-line-number="1707"></td>
        <td id="LC1707" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1708" class="blob-num js-line-number" data-line-number="1708"></td>
        <td id="LC1708" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1709" class="blob-num js-line-number" data-line-number="1709"></td>
        <td id="LC1709" class="blob-code js-file-line">        }  <span class="pl-c">//for iAzim</span></td>
      </tr>
      <tr>
        <td id="L1710" class="blob-num js-line-number" data-line-number="1710"></td>
        <td id="LC1710" class="blob-code js-file-line">    } <span class="pl-c">//for iRang</span></td>
      </tr>
      <tr>
        <td id="L1711" class="blob-num js-line-number" data-line-number="1711"></td>
        <td id="LC1711" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1712" class="blob-num js-line-number" data-line-number="1712"></td>
        <td id="LC1712" class="blob-code js-file-line">    <span class="pl-k">return</span> nPointsWritten;</td>
      </tr>
      <tr>
        <td id="L1713" class="blob-num js-line-number" data-line-number="1713"></td>
        <td id="LC1713" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1714" class="blob-num js-line-number" data-line-number="1714"></td>
        <td id="LC1714" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1715" class="blob-num js-line-number" data-line-number="1715"></td>
        <td id="LC1715" class="blob-code js-file-line">} <span class="pl-c">// getListOfSelectedGates</span></td>
      </tr>
      <tr>
        <td id="L1716" class="blob-num js-line-number" data-line-number="1716"></td>
        <td id="LC1716" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1717" class="blob-num js-line-number" data-line-number="1717"></td>
        <td id="LC1717" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1718" class="blob-num js-line-number" data-line-number="1718"></td>
        <td id="LC1718" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1719" class="blob-num js-line-number" data-line-number="1719"></td>
        <td id="LC1719" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1720" class="blob-num js-line-number" data-line-number="1720"></td>
        <td id="LC1720" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">hasAzimuthGap</span>(<span class="pl-k">const</span> <span class="pl-k">float</span>* points, <span class="pl-k">const</span> <span class="pl-k">int</span> nPoints) {</td>
      </tr>
      <tr>
        <td id="L1721" class="blob-num js-line-number" data-line-number="1721"></td>
        <td id="LC1721" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1722" class="blob-num js-line-number" data-line-number="1722"></td>
        <td id="LC1722" class="blob-code js-file-line">    <span class="pl-k">int</span> hasGap;</td>
      </tr>
      <tr>
        <td id="L1723" class="blob-num js-line-number" data-line-number="1723"></td>
        <td id="LC1723" class="blob-code js-file-line">    <span class="pl-k">int</span> nObs[nBinsGap];</td>
      </tr>
      <tr>
        <td id="L1724" class="blob-num js-line-number" data-line-number="1724"></td>
        <td id="LC1724" class="blob-code js-file-line">    <span class="pl-k">int</span> iPoint;</td>
      </tr>
      <tr>
        <td id="L1725" class="blob-num js-line-number" data-line-number="1725"></td>
        <td id="LC1725" class="blob-code js-file-line">    <span class="pl-k">int</span> iBinGap;</td>
      </tr>
      <tr>
        <td id="L1726" class="blob-num js-line-number" data-line-number="1726"></td>
        <td id="LC1726" class="blob-code js-file-line">    <span class="pl-k">int</span> iBinGapNext;</td>
      </tr>
      <tr>
        <td id="L1727" class="blob-num js-line-number" data-line-number="1727"></td>
        <td id="LC1727" class="blob-code js-file-line">    <span class="pl-k">float</span> azimuth;</td>
      </tr>
      <tr>
        <td id="L1728" class="blob-num js-line-number" data-line-number="1728"></td>
        <td id="LC1728" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1729" class="blob-num js-line-number" data-line-number="1729"></td>
        <td id="LC1729" class="blob-code js-file-line">    hasGap = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1730" class="blob-num js-line-number" data-line-number="1730"></td>
        <td id="LC1730" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1731" class="blob-num js-line-number" data-line-number="1731"></td>
        <td id="LC1731" class="blob-code js-file-line">    <span class="pl-c">// Initialize histogram</span></td>
      </tr>
      <tr>
        <td id="L1732" class="blob-num js-line-number" data-line-number="1732"></td>
        <td id="LC1732" class="blob-code js-file-line">    <span class="pl-k">for</span> (iBinGap = <span class="pl-c1">0</span>; iBinGap &lt; nBinsGap; iBinGap++) {</td>
      </tr>
      <tr>
        <td id="L1733" class="blob-num js-line-number" data-line-number="1733"></td>
        <td id="LC1733" class="blob-code js-file-line">        nObs[iBinGap] = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L1734" class="blob-num js-line-number" data-line-number="1734"></td>
        <td id="LC1734" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1735" class="blob-num js-line-number" data-line-number="1735"></td>
        <td id="LC1735" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1736" class="blob-num js-line-number" data-line-number="1736"></td>
        <td id="LC1736" class="blob-code js-file-line">    <span class="pl-c">// Collect histogram data</span></td>
      </tr>
      <tr>
        <td id="L1737" class="blob-num js-line-number" data-line-number="1737"></td>
        <td id="LC1737" class="blob-code js-file-line">    <span class="pl-k">for</span> (iPoint = <span class="pl-c1">0</span>; iPoint &lt; nPoints; iPoint++) {</td>
      </tr>
      <tr>
        <td id="L1738" class="blob-num js-line-number" data-line-number="1738"></td>
        <td id="LC1738" class="blob-code js-file-line">        azimuth = points[iPoint*nDims];</td>
      </tr>
      <tr>
        <td id="L1739" class="blob-num js-line-number" data-line-number="1739"></td>
        <td id="LC1739" class="blob-code js-file-line">        iBinGap = ((<span class="pl-k">int</span>) <span class="pl-c1">floor</span>((azimuth / <span class="pl-c1">360.0</span>) * nBinsGap)) % nBinsGap;</td>
      </tr>
      <tr>
        <td id="L1740" class="blob-num js-line-number" data-line-number="1740"></td>
        <td id="LC1740" class="blob-code js-file-line">        nObs[iBinGap]++;</td>
      </tr>
      <tr>
        <td id="L1741" class="blob-num js-line-number" data-line-number="1741"></td>
        <td id="LC1741" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1742" class="blob-num js-line-number" data-line-number="1742"></td>
        <td id="LC1742" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1743" class="blob-num js-line-number" data-line-number="1743"></td>
        <td id="LC1743" class="blob-code js-file-line">    <span class="pl-c">// Detect adjacent bins in which the number of azimuth observations </span></td>
      </tr>
      <tr>
        <td id="L1744" class="blob-num js-line-number" data-line-number="1744"></td>
        <td id="LC1744" class="blob-code js-file-line">    <span class="pl-c">// is less than the minimum required number</span></td>
      </tr>
      <tr>
        <td id="L1745" class="blob-num js-line-number" data-line-number="1745"></td>
        <td id="LC1745" class="blob-code js-file-line">    <span class="pl-k">for</span> (iBinGap = <span class="pl-c1">0</span>; iBinGap &lt; nBinsGap; iBinGap++) {</td>
      </tr>
      <tr>
        <td id="L1746" class="blob-num js-line-number" data-line-number="1746"></td>
        <td id="LC1746" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1747" class="blob-num js-line-number" data-line-number="1747"></td>
        <td id="LC1747" class="blob-code js-file-line">        iBinGapNext = (iBinGap + <span class="pl-c1">1</span>) % nBinsGap;</td>
      </tr>
      <tr>
        <td id="L1748" class="blob-num js-line-number" data-line-number="1748"></td>
        <td id="LC1748" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1749" class="blob-num js-line-number" data-line-number="1749"></td>
        <td id="LC1749" class="blob-code js-file-line">        <span class="pl-k">if</span> (nObs[iBinGap] &lt; nObsGapMin &amp;&amp; nObs[iBinGapNext] &lt; nObsGapMin) {</td>
      </tr>
      <tr>
        <td id="L1750" class="blob-num js-line-number" data-line-number="1750"></td>
        <td id="LC1750" class="blob-code js-file-line">            hasGap = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1751" class="blob-num js-line-number" data-line-number="1751"></td>
        <td id="LC1751" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1752" class="blob-num js-line-number" data-line-number="1752"></td>
        <td id="LC1752" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1753" class="blob-num js-line-number" data-line-number="1753"></td>
        <td id="LC1753" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1754" class="blob-num js-line-number" data-line-number="1754"></td>
        <td id="LC1754" class="blob-code js-file-line">    <span class="pl-k">return</span> hasGap;</td>
      </tr>
      <tr>
        <td id="L1755" class="blob-num js-line-number" data-line-number="1755"></td>
        <td id="LC1755" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1756" class="blob-num js-line-number" data-line-number="1756"></td>
        <td id="LC1756" class="blob-code js-file-line">} <span class="pl-c">// hasAzimuthGap</span></td>
      </tr>
      <tr>
        <td id="L1757" class="blob-num js-line-number" data-line-number="1757"></td>
        <td id="LC1757" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1758" class="blob-num js-line-number" data-line-number="1758"></td>
        <td id="LC1758" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1759" class="blob-num js-line-number" data-line-number="1759"></td>
        <td id="LC1759" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1760" class="blob-num js-line-number" data-line-number="1760"></td>
        <td id="LC1760" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1761" class="blob-num js-line-number" data-line-number="1761"></td>
        <td id="LC1761" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1762" class="blob-num js-line-number" data-line-number="1762"></td>
        <td id="LC1762" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1763" class="blob-num js-line-number" data-line-number="1763"></td>
        <td id="LC1763" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">includeGate</span>(<span class="pl-k">const</span> <span class="pl-k">int</span> iProfileType, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">int</span> gateCode) {</td>
      </tr>
      <tr>
        <td id="L1764" class="blob-num js-line-number" data-line-number="1764"></td>
        <td id="LC1764" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1765" class="blob-num js-line-number" data-line-number="1765"></td>
        <td id="LC1765" class="blob-code js-file-line">    <span class="pl-k">int</span> doInclude = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1766" class="blob-num js-line-number" data-line-number="1766"></td>
        <td id="LC1766" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1767" class="blob-num js-line-number" data-line-number="1767"></td>
        <td id="LC1767" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionStaticClutter) {</td>
      </tr>
      <tr>
        <td id="L1768" class="blob-num js-line-number" data-line-number="1768"></td>
        <td id="LC1768" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1769" class="blob-num js-line-number" data-line-number="1769"></td>
        <td id="LC1769" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 0 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1770" class="blob-num js-line-number" data-line-number="1770"></td>
        <td id="LC1770" class="blob-code js-file-line">        <span class="pl-c">// this gate is true in the static clutter map (which we don&#39;t have yet TODO)</span></td>
      </tr>
      <tr>
        <td id="L1771" class="blob-num js-line-number" data-line-number="1771"></td>
        <td id="LC1771" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1772" class="blob-num js-line-number" data-line-number="1772"></td>
        <td id="LC1772" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1773" class="blob-num js-line-number" data-line-number="1773"></td>
        <td id="LC1773" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1774" class="blob-num js-line-number" data-line-number="1774"></td>
        <td id="LC1774" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1775" class="blob-num js-line-number" data-line-number="1775"></td>
        <td id="LC1775" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1776" class="blob-num js-line-number" data-line-number="1776"></td>
        <td id="LC1776" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1777" class="blob-num js-line-number" data-line-number="1777"></td>
        <td id="LC1777" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1778" class="blob-num js-line-number" data-line-number="1778"></td>
        <td id="LC1778" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1779" class="blob-num js-line-number" data-line-number="1779"></td>
        <td id="LC1779" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1780" class="blob-num js-line-number" data-line-number="1780"></td>
        <td id="LC1780" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1781" class="blob-num js-line-number" data-line-number="1781"></td>
        <td id="LC1781" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1782" class="blob-num js-line-number" data-line-number="1782"></td>
        <td id="LC1782" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1783" class="blob-num js-line-number" data-line-number="1783"></td>
        <td id="LC1783" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1784" class="blob-num js-line-number" data-line-number="1784"></td>
        <td id="LC1784" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1785" class="blob-num js-line-number" data-line-number="1785"></td>
        <td id="LC1785" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1786" class="blob-num js-line-number" data-line-number="1786"></td>
        <td id="LC1786" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1787" class="blob-num js-line-number" data-line-number="1787"></td>
        <td id="LC1787" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionDynamicClutter) {</td>
      </tr>
      <tr>
        <td id="L1788" class="blob-num js-line-number" data-line-number="1788"></td>
        <td id="LC1788" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1789" class="blob-num js-line-number" data-line-number="1789"></td>
        <td id="LC1789" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 1 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1790" class="blob-num js-line-number" data-line-number="1790"></td>
        <td id="LC1790" class="blob-code js-file-line">        <span class="pl-c">// this gate is part of the cluttermap (without fringe)</span></td>
      </tr>
      <tr>
        <td id="L1791" class="blob-num js-line-number" data-line-number="1791"></td>
        <td id="LC1791" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1792" class="blob-num js-line-number" data-line-number="1792"></td>
        <td id="LC1792" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1793" class="blob-num js-line-number" data-line-number="1793"></td>
        <td id="LC1793" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1794" class="blob-num js-line-number" data-line-number="1794"></td>
        <td id="LC1794" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1795" class="blob-num js-line-number" data-line-number="1795"></td>
        <td id="LC1795" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1796" class="blob-num js-line-number" data-line-number="1796"></td>
        <td id="LC1796" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1797" class="blob-num js-line-number" data-line-number="1797"></td>
        <td id="LC1797" class="blob-code js-file-line">                doInclude = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1798" class="blob-num js-line-number" data-line-number="1798"></td>
        <td id="LC1798" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1799" class="blob-num js-line-number" data-line-number="1799"></td>
        <td id="LC1799" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1800" class="blob-num js-line-number" data-line-number="1800"></td>
        <td id="LC1800" class="blob-code js-file-line">                doInclude = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1801" class="blob-num js-line-number" data-line-number="1801"></td>
        <td id="LC1801" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1802" class="blob-num js-line-number" data-line-number="1802"></td>
        <td id="LC1802" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1803" class="blob-num js-line-number" data-line-number="1803"></td>
        <td id="LC1803" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1804" class="blob-num js-line-number" data-line-number="1804"></td>
        <td id="LC1804" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1805" class="blob-num js-line-number" data-line-number="1805"></td>
        <td id="LC1805" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1806" class="blob-num js-line-number" data-line-number="1806"></td>
        <td id="LC1806" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1807" class="blob-num js-line-number" data-line-number="1807"></td>
        <td id="LC1807" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionDynamicClutterFringe) {</td>
      </tr>
      <tr>
        <td id="L1808" class="blob-num js-line-number" data-line-number="1808"></td>
        <td id="LC1808" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1809" class="blob-num js-line-number" data-line-number="1809"></td>
        <td id="LC1809" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 2 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1810" class="blob-num js-line-number" data-line-number="1810"></td>
        <td id="LC1810" class="blob-code js-file-line">        <span class="pl-c">// this gate is part of the fringe of the cluttermap</span></td>
      </tr>
      <tr>
        <td id="L1811" class="blob-num js-line-number" data-line-number="1811"></td>
        <td id="LC1811" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1812" class="blob-num js-line-number" data-line-number="1812"></td>
        <td id="LC1812" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1813" class="blob-num js-line-number" data-line-number="1813"></td>
        <td id="LC1813" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1814" class="blob-num js-line-number" data-line-number="1814"></td>
        <td id="LC1814" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1815" class="blob-num js-line-number" data-line-number="1815"></td>
        <td id="LC1815" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1816" class="blob-num js-line-number" data-line-number="1816"></td>
        <td id="LC1816" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1817" class="blob-num js-line-number" data-line-number="1817"></td>
        <td id="LC1817" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1818" class="blob-num js-line-number" data-line-number="1818"></td>
        <td id="LC1818" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1819" class="blob-num js-line-number" data-line-number="1819"></td>
        <td id="LC1819" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1820" class="blob-num js-line-number" data-line-number="1820"></td>
        <td id="LC1820" class="blob-code js-file-line">                doInclude = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1821" class="blob-num js-line-number" data-line-number="1821"></td>
        <td id="LC1821" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1822" class="blob-num js-line-number" data-line-number="1822"></td>
        <td id="LC1822" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1823" class="blob-num js-line-number" data-line-number="1823"></td>
        <td id="LC1823" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1824" class="blob-num js-line-number" data-line-number="1824"></td>
        <td id="LC1824" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1825" class="blob-num js-line-number" data-line-number="1825"></td>
        <td id="LC1825" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1826" class="blob-num js-line-number" data-line-number="1826"></td>
        <td id="LC1826" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1827" class="blob-num js-line-number" data-line-number="1827"></td>
        <td id="LC1827" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionVradMissing) {</td>
      </tr>
      <tr>
        <td id="L1828" class="blob-num js-line-number" data-line-number="1828"></td>
        <td id="LC1828" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1829" class="blob-num js-line-number" data-line-number="1829"></td>
        <td id="LC1829" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 3 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1830" class="blob-num js-line-number" data-line-number="1830"></td>
        <td id="LC1830" class="blob-code js-file-line">        <span class="pl-c">// this gate has reflectivity data but no corresponding radial velocity data</span></td>
      </tr>
      <tr>
        <td id="L1831" class="blob-num js-line-number" data-line-number="1831"></td>
        <td id="LC1831" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1832" class="blob-num js-line-number" data-line-number="1832"></td>
        <td id="LC1832" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1833" class="blob-num js-line-number" data-line-number="1833"></td>
        <td id="LC1833" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1834" class="blob-num js-line-number" data-line-number="1834"></td>
        <td id="LC1834" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1835" class="blob-num js-line-number" data-line-number="1835"></td>
        <td id="LC1835" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1836" class="blob-num js-line-number" data-line-number="1836"></td>
        <td id="LC1836" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1837" class="blob-num js-line-number" data-line-number="1837"></td>
        <td id="LC1837" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1838" class="blob-num js-line-number" data-line-number="1838"></td>
        <td id="LC1838" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1839" class="blob-num js-line-number" data-line-number="1839"></td>
        <td id="LC1839" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1840" class="blob-num js-line-number" data-line-number="1840"></td>
        <td id="LC1840" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1841" class="blob-num js-line-number" data-line-number="1841"></td>
        <td id="LC1841" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1842" class="blob-num js-line-number" data-line-number="1842"></td>
        <td id="LC1842" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1843" class="blob-num js-line-number" data-line-number="1843"></td>
        <td id="LC1843" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1844" class="blob-num js-line-number" data-line-number="1844"></td>
        <td id="LC1844" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1845" class="blob-num js-line-number" data-line-number="1845"></td>
        <td id="LC1845" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1846" class="blob-num js-line-number" data-line-number="1846"></td>
        <td id="LC1846" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1847" class="blob-num js-line-number" data-line-number="1847"></td>
        <td id="LC1847" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionDbzTooHighForBirds) {</td>
      </tr>
      <tr>
        <td id="L1848" class="blob-num js-line-number" data-line-number="1848"></td>
        <td id="LC1848" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1849" class="blob-num js-line-number" data-line-number="1849"></td>
        <td id="LC1849" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 4 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1850" class="blob-num js-line-number" data-line-number="1850"></td>
        <td id="LC1850" class="blob-code js-file-line">        <span class="pl-c">// this gate&#39;s dbz value is too high to be due to birds, it must be </span></td>
      </tr>
      <tr>
        <td id="L1851" class="blob-num js-line-number" data-line-number="1851"></td>
        <td id="LC1851" class="blob-code js-file-line">        <span class="pl-c">// caused by something else</span></td>
      </tr>
      <tr>
        <td id="L1852" class="blob-num js-line-number" data-line-number="1852"></td>
        <td id="LC1852" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1853" class="blob-num js-line-number" data-line-number="1853"></td>
        <td id="LC1853" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1854" class="blob-num js-line-number" data-line-number="1854"></td>
        <td id="LC1854" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1855" class="blob-num js-line-number" data-line-number="1855"></td>
        <td id="LC1855" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1856" class="blob-num js-line-number" data-line-number="1856"></td>
        <td id="LC1856" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1857" class="blob-num js-line-number" data-line-number="1857"></td>
        <td id="LC1857" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1858" class="blob-num js-line-number" data-line-number="1858"></td>
        <td id="LC1858" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1859" class="blob-num js-line-number" data-line-number="1859"></td>
        <td id="LC1859" class="blob-code js-file-line">                doInclude = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1860" class="blob-num js-line-number" data-line-number="1860"></td>
        <td id="LC1860" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1861" class="blob-num js-line-number" data-line-number="1861"></td>
        <td id="LC1861" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1862" class="blob-num js-line-number" data-line-number="1862"></td>
        <td id="LC1862" class="blob-code js-file-line">                doInclude = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1863" class="blob-num js-line-number" data-line-number="1863"></td>
        <td id="LC1863" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1864" class="blob-num js-line-number" data-line-number="1864"></td>
        <td id="LC1864" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1865" class="blob-num js-line-number" data-line-number="1865"></td>
        <td id="LC1865" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1866" class="blob-num js-line-number" data-line-number="1866"></td>
        <td id="LC1866" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1867" class="blob-num js-line-number" data-line-number="1867"></td>
        <td id="LC1867" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1868" class="blob-num js-line-number" data-line-number="1868"></td>
        <td id="LC1868" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1869" class="blob-num js-line-number" data-line-number="1869"></td>
        <td id="LC1869" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionVradTooLow) {</td>
      </tr>
      <tr>
        <td id="L1870" class="blob-num js-line-number" data-line-number="1870"></td>
        <td id="LC1870" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1871" class="blob-num js-line-number" data-line-number="1871"></td>
        <td id="LC1871" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 5 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1872" class="blob-num js-line-number" data-line-number="1872"></td>
        <td id="LC1872" class="blob-code js-file-line">        <span class="pl-c">// this gate&#39;s radial velocity is too low to be due to actual scatterers; likely just noise</span></td>
      </tr>
      <tr>
        <td id="L1873" class="blob-num js-line-number" data-line-number="1873"></td>
        <td id="LC1873" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1874" class="blob-num js-line-number" data-line-number="1874"></td>
        <td id="LC1874" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1875" class="blob-num js-line-number" data-line-number="1875"></td>
        <td id="LC1875" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1876" class="blob-num js-line-number" data-line-number="1876"></td>
        <td id="LC1876" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1877" class="blob-num js-line-number" data-line-number="1877"></td>
        <td id="LC1877" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1878" class="blob-num js-line-number" data-line-number="1878"></td>
        <td id="LC1878" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1879" class="blob-num js-line-number" data-line-number="1879"></td>
        <td id="LC1879" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1880" class="blob-num js-line-number" data-line-number="1880"></td>
        <td id="LC1880" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1881" class="blob-num js-line-number" data-line-number="1881"></td>
        <td id="LC1881" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1882" class="blob-num js-line-number" data-line-number="1882"></td>
        <td id="LC1882" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1883" class="blob-num js-line-number" data-line-number="1883"></td>
        <td id="LC1883" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1884" class="blob-num js-line-number" data-line-number="1884"></td>
        <td id="LC1884" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1885" class="blob-num js-line-number" data-line-number="1885"></td>
        <td id="LC1885" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1886" class="blob-num js-line-number" data-line-number="1886"></td>
        <td id="LC1886" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1887" class="blob-num js-line-number" data-line-number="1887"></td>
        <td id="LC1887" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1888" class="blob-num js-line-number" data-line-number="1888"></td>
        <td id="LC1888" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1889" class="blob-num js-line-number" data-line-number="1889"></td>
        <td id="LC1889" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionVDifMax) {</td>
      </tr>
      <tr>
        <td id="L1890" class="blob-num js-line-number" data-line-number="1890"></td>
        <td id="LC1890" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1891" class="blob-num js-line-number" data-line-number="1891"></td>
        <td id="LC1891" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 6 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1892" class="blob-num js-line-number" data-line-number="1892"></td>
        <td id="LC1892" class="blob-code js-file-line">        <span class="pl-c">// after the first svdfit, this gate&#39;s fitted vRad was more than </span></td>
      </tr>
      <tr>
        <td id="L1893" class="blob-num js-line-number" data-line-number="1893"></td>
        <td id="LC1893" class="blob-code js-file-line">        <span class="pl-c">// VDIFMAX away from the observed vRad for that gate. It is therefore</span></td>
      </tr>
      <tr>
        <td id="L1894" class="blob-num js-line-number" data-line-number="1894"></td>
        <td id="LC1894" class="blob-code js-file-line">        <span class="pl-c">// considered an outlier</span></td>
      </tr>
      <tr>
        <td id="L1895" class="blob-num js-line-number" data-line-number="1895"></td>
        <td id="LC1895" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1896" class="blob-num js-line-number" data-line-number="1896"></td>
        <td id="LC1896" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1897" class="blob-num js-line-number" data-line-number="1897"></td>
        <td id="LC1897" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1898" class="blob-num js-line-number" data-line-number="1898"></td>
        <td id="LC1898" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1899" class="blob-num js-line-number" data-line-number="1899"></td>
        <td id="LC1899" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1900" class="blob-num js-line-number" data-line-number="1900"></td>
        <td id="LC1900" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1901" class="blob-num js-line-number" data-line-number="1901"></td>
        <td id="LC1901" class="blob-code js-file-line">                <span class="pl-c">// FIXME not sure if this should be FALSE for iProfileType==2</span></td>
      </tr>
      <tr>
        <td id="L1902" class="blob-num js-line-number" data-line-number="1902"></td>
        <td id="LC1902" class="blob-code js-file-line">                doInclude = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1903" class="blob-num js-line-number" data-line-number="1903"></td>
        <td id="LC1903" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1904" class="blob-num js-line-number" data-line-number="1904"></td>
        <td id="LC1904" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1905" class="blob-num js-line-number" data-line-number="1905"></td>
        <td id="LC1905" class="blob-code js-file-line">                <span class="pl-c">// FIXME not sure if this should be FALSE for iProfileType==3</span></td>
      </tr>
      <tr>
        <td id="L1906" class="blob-num js-line-number" data-line-number="1906"></td>
        <td id="LC1906" class="blob-code js-file-line">                doInclude = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L1907" class="blob-num js-line-number" data-line-number="1907"></td>
        <td id="LC1907" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1908" class="blob-num js-line-number" data-line-number="1908"></td>
        <td id="LC1908" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1909" class="blob-num js-line-number" data-line-number="1909"></td>
        <td id="LC1909" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1910" class="blob-num js-line-number" data-line-number="1910"></td>
        <td id="LC1910" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1911" class="blob-num js-line-number" data-line-number="1911"></td>
        <td id="LC1911" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1912" class="blob-num js-line-number" data-line-number="1912"></td>
        <td id="LC1912" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1913" class="blob-num js-line-number" data-line-number="1913"></td>
        <td id="LC1913" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1914" class="blob-num js-line-number" data-line-number="1914"></td>
        <td id="LC1914" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1915" class="blob-num js-line-number" data-line-number="1915"></td>
        <td id="LC1915" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionAzimTooLow) {</td>
      </tr>
      <tr>
        <td id="L1916" class="blob-num js-line-number" data-line-number="1916"></td>
        <td id="LC1916" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1917" class="blob-num js-line-number" data-line-number="1917"></td>
        <td id="LC1917" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 7 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1918" class="blob-num js-line-number" data-line-number="1918"></td>
        <td id="LC1918" class="blob-code js-file-line">        <span class="pl-c">// the user can specify to exclude gates based on their azimuth;</span></td>
      </tr>
      <tr>
        <td id="L1919" class="blob-num js-line-number" data-line-number="1919"></td>
        <td id="LC1919" class="blob-code js-file-line">        <span class="pl-c">// this clause is for gates that have too low azimuth</span></td>
      </tr>
      <tr>
        <td id="L1920" class="blob-num js-line-number" data-line-number="1920"></td>
        <td id="LC1920" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1921" class="blob-num js-line-number" data-line-number="1921"></td>
        <td id="LC1921" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1922" class="blob-num js-line-number" data-line-number="1922"></td>
        <td id="LC1922" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1923" class="blob-num js-line-number" data-line-number="1923"></td>
        <td id="LC1923" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1924" class="blob-num js-line-number" data-line-number="1924"></td>
        <td id="LC1924" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1925" class="blob-num js-line-number" data-line-number="1925"></td>
        <td id="LC1925" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1926" class="blob-num js-line-number" data-line-number="1926"></td>
        <td id="LC1926" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1927" class="blob-num js-line-number" data-line-number="1927"></td>
        <td id="LC1927" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1928" class="blob-num js-line-number" data-line-number="1928"></td>
        <td id="LC1928" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1929" class="blob-num js-line-number" data-line-number="1929"></td>
        <td id="LC1929" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1930" class="blob-num js-line-number" data-line-number="1930"></td>
        <td id="LC1930" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1931" class="blob-num js-line-number" data-line-number="1931"></td>
        <td id="LC1931" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1932" class="blob-num js-line-number" data-line-number="1932"></td>
        <td id="LC1932" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1933" class="blob-num js-line-number" data-line-number="1933"></td>
        <td id="LC1933" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1934" class="blob-num js-line-number" data-line-number="1934"></td>
        <td id="LC1934" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1935" class="blob-num js-line-number" data-line-number="1935"></td>
        <td id="LC1935" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1936" class="blob-num js-line-number" data-line-number="1936"></td>
        <td id="LC1936" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1937" class="blob-num js-line-number" data-line-number="1937"></td>
        <td id="LC1937" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &amp; <span class="pl-c1">1</span>&lt;&lt;flagPositionAzimTooHigh) {</td>
      </tr>
      <tr>
        <td id="L1938" class="blob-num js-line-number" data-line-number="1938"></td>
        <td id="LC1938" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1939" class="blob-num js-line-number" data-line-number="1939"></td>
        <td id="LC1939" class="blob-code js-file-line">        <span class="pl-c">// i.e. flag 8 in gateCode is true</span></td>
      </tr>
      <tr>
        <td id="L1940" class="blob-num js-line-number" data-line-number="1940"></td>
        <td id="LC1940" class="blob-code js-file-line">        <span class="pl-c">// the user can specify to exclude gates based on their azimuth;</span></td>
      </tr>
      <tr>
        <td id="L1941" class="blob-num js-line-number" data-line-number="1941"></td>
        <td id="LC1941" class="blob-code js-file-line">        <span class="pl-c">// this clause is for gates that have too high azimuth</span></td>
      </tr>
      <tr>
        <td id="L1942" class="blob-num js-line-number" data-line-number="1942"></td>
        <td id="LC1942" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L1943" class="blob-num js-line-number" data-line-number="1943"></td>
        <td id="LC1943" class="blob-code js-file-line">        <span class="pl-k">switch</span> (iProfileType) {</td>
      </tr>
      <tr>
        <td id="L1944" class="blob-num js-line-number" data-line-number="1944"></td>
        <td id="LC1944" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">1</span> : </td>
      </tr>
      <tr>
        <td id="L1945" class="blob-num js-line-number" data-line-number="1945"></td>
        <td id="LC1945" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1946" class="blob-num js-line-number" data-line-number="1946"></td>
        <td id="LC1946" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1947" class="blob-num js-line-number" data-line-number="1947"></td>
        <td id="LC1947" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">2</span> : </td>
      </tr>
      <tr>
        <td id="L1948" class="blob-num js-line-number" data-line-number="1948"></td>
        <td id="LC1948" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1949" class="blob-num js-line-number" data-line-number="1949"></td>
        <td id="LC1949" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1950" class="blob-num js-line-number" data-line-number="1950"></td>
        <td id="LC1950" class="blob-code js-file-line">            <span class="pl-k">case</span> <span class="pl-c1">3</span> : </td>
      </tr>
      <tr>
        <td id="L1951" class="blob-num js-line-number" data-line-number="1951"></td>
        <td id="LC1951" class="blob-code js-file-line">                doInclude = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L1952" class="blob-num js-line-number" data-line-number="1952"></td>
        <td id="LC1952" class="blob-code js-file-line">                <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L1953" class="blob-num js-line-number" data-line-number="1953"></td>
        <td id="LC1953" class="blob-code js-file-line">            <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L1954" class="blob-num js-line-number" data-line-number="1954"></td>
        <td id="LC1954" class="blob-code js-file-line">                <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>Something went wrong; behavior not implemented for given iProfileType.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L1955" class="blob-num js-line-number" data-line-number="1955"></td>
        <td id="LC1955" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L1956" class="blob-num js-line-number" data-line-number="1956"></td>
        <td id="LC1956" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1957" class="blob-num js-line-number" data-line-number="1957"></td>
        <td id="LC1957" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1958" class="blob-num js-line-number" data-line-number="1958"></td>
        <td id="LC1958" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1959" class="blob-num js-line-number" data-line-number="1959"></td>
        <td id="LC1959" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1960" class="blob-num js-line-number" data-line-number="1960"></td>
        <td id="LC1960" class="blob-code js-file-line">    <span class="pl-k">return</span> doInclude;</td>
      </tr>
      <tr>
        <td id="L1961" class="blob-num js-line-number" data-line-number="1961"></td>
        <td id="LC1961" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1962" class="blob-num js-line-number" data-line-number="1962"></td>
        <td id="LC1962" class="blob-code js-file-line">} <span class="pl-c">// includeGate</span></td>
      </tr>
      <tr>
        <td id="L1963" class="blob-num js-line-number" data-line-number="1963"></td>
        <td id="LC1963" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1964" class="blob-num js-line-number" data-line-number="1964"></td>
        <td id="LC1964" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1965" class="blob-num js-line-number" data-line-number="1965"></td>
        <td id="LC1965" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1966" class="blob-num js-line-number" data-line-number="1966"></td>
        <td id="LC1966" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1967" class="blob-num js-line-number" data-line-number="1967"></td>
        <td id="LC1967" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1968" class="blob-num js-line-number" data-line-number="1968"></td>
        <td id="LC1968" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">readUserConfigOptions</span>(<span class="pl-k">void</span>) {</td>
      </tr>
      <tr>
        <td id="L1969" class="blob-num js-line-number" data-line-number="1969"></td>
        <td id="LC1969" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1970" class="blob-num js-line-number" data-line-number="1970"></td>
        <td id="LC1970" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L1971" class="blob-num js-line-number" data-line-number="1971"></td>
        <td id="LC1971" class="blob-code js-file-line">    <span class="pl-c1">cfg_opt_t</span> opts[] = {</td>
      </tr>
      <tr>
        <td id="L1972" class="blob-num js-line-number" data-line-number="1972"></td>
        <td id="LC1972" class="blob-code js-file-line">        <span class="pl-c1">CFG_FLOAT</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>HLAYER<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">200</span>.<span class="pl-c1">0f</span>, CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1973" class="blob-num js-line-number" data-line-number="1973"></td>
        <td id="LC1973" class="blob-code js-file-line">        <span class="pl-c1">CFG_INT</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>NLAYER<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">30</span>, CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1974" class="blob-num js-line-number" data-line-number="1974"></td>
        <td id="LC1974" class="blob-code js-file-line">        <span class="pl-c1">CFG_FLOAT</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>RANGEMIN<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">5000</span>.<span class="pl-c1">0f</span>, CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1975" class="blob-num js-line-number" data-line-number="1975"></td>
        <td id="LC1975" class="blob-code js-file-line">        <span class="pl-c1">CFG_FLOAT</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>RANGEMAX<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">25000</span>.<span class="pl-c1">0f</span>, CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1976" class="blob-num js-line-number" data-line-number="1976"></td>
        <td id="LC1976" class="blob-code js-file-line">        <span class="pl-c1">CFG_FLOAT</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>AZIMMIN<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">0</span>.<span class="pl-c1">0f</span>, CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1977" class="blob-num js-line-number" data-line-number="1977"></td>
        <td id="LC1977" class="blob-code js-file-line">        <span class="pl-c1">CFG_FLOAT</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>AZIMMAX<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">360</span>.<span class="pl-c1">0f</span>, CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1978" class="blob-num js-line-number" data-line-number="1978"></td>
        <td id="LC1978" class="blob-code js-file-line">        <span class="pl-c1">CFG_FLOAT</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>RADAR_WAVELENGTH_CM<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">5</span>.<span class="pl-c1">3f</span>, CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1979" class="blob-num js-line-number" data-line-number="1979"></td>
        <td id="LC1979" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>USE_STATIC_CLUTTER_DATA<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">FALSE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1980" class="blob-num js-line-number" data-line-number="1980"></td>
        <td id="LC1980" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>VERBOSE_OUTPUT_REQUIRED<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">FALSE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1981" class="blob-num js-line-number" data-line-number="1981"></td>
        <td id="LC1981" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_DBZ<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1982" class="blob-num js-line-number" data-line-number="1982"></td>
        <td id="LC1982" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_VRAD<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1983" class="blob-num js-line-number" data-line-number="1983"></td>
        <td id="LC1983" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_CELL<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1984" class="blob-num js-line-number" data-line-number="1984"></td>
        <td id="LC1984" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_CELL_PROP<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1985" class="blob-num js-line-number" data-line-number="1985"></td>
        <td id="LC1985" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_TEXTURE<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1986" class="blob-num js-line-number" data-line-number="1986"></td>
        <td id="LC1986" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_CLUT<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1987" class="blob-num js-line-number" data-line-number="1987"></td>
        <td id="LC1987" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_OPTIONS<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1988" class="blob-num js-line-number" data-line-number="1988"></td>
        <td id="LC1988" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>FIT_VRAD<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1989" class="blob-num js-line-number" data-line-number="1989"></td>
        <td id="LC1989" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_PROFILE<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">TRUE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1990" class="blob-num js-line-number" data-line-number="1990"></td>
        <td id="LC1990" class="blob-code js-file-line">        <span class="pl-c1">CFG_BOOL</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_POINTS_ARRAY<span class="pl-pds">&quot;</span></span>,<span class="pl-c1">FALSE</span>,CFGF_NONE),</td>
      </tr>
      <tr>
        <td id="L1991" class="blob-num js-line-number" data-line-number="1991"></td>
        <td id="LC1991" class="blob-code js-file-line">        <span class="pl-c1">CFG_END</span>()</td>
      </tr>
      <tr>
        <td id="L1992" class="blob-num js-line-number" data-line-number="1992"></td>
        <td id="LC1992" class="blob-code js-file-line">    };</td>
      </tr>
      <tr>
        <td id="L1993" class="blob-num js-line-number" data-line-number="1993"></td>
        <td id="LC1993" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1994" class="blob-num js-line-number" data-line-number="1994"></td>
        <td id="LC1994" class="blob-code js-file-line">    cfg = <span class="pl-c1">cfg_init</span>(opts, CFGF_NONE);</td>
      </tr>
      <tr>
        <td id="L1995" class="blob-num js-line-number" data-line-number="1995"></td>
        <td id="LC1995" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L1996" class="blob-num js-line-number" data-line-number="1996"></td>
        <td id="LC1996" class="blob-code js-file-line">    <span class="pl-k">if</span> (<span class="pl-c1">cfg_parse</span>(cfg, <span class="pl-s"><span class="pl-pds">&quot;</span>options.conf<span class="pl-pds">&quot;</span></span>) == CFG_PARSE_ERROR) {</td>
      </tr>
      <tr>
        <td id="L1997" class="blob-num js-line-number" data-line-number="1997"></td>
        <td id="LC1997" class="blob-code js-file-line">        <span class="pl-k">return</span> <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L1998" class="blob-num js-line-number" data-line-number="1998"></td>
        <td id="LC1998" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L1999" class="blob-num js-line-number" data-line-number="1999"></td>
        <td id="LC1999" class="blob-code js-file-line">   </td>
      </tr>
      <tr>
        <td id="L2000" class="blob-num js-line-number" data-line-number="2000"></td>
        <td id="LC2000" class="blob-code js-file-line">    <span class="pl-k">return</span> <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2001" class="blob-num js-line-number" data-line-number="2001"></td>
        <td id="LC2001" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2002" class="blob-num js-line-number" data-line-number="2002"></td>
        <td id="LC2002" class="blob-code js-file-line">} <span class="pl-c">// readUserConfigOptions</span></td>
      </tr>
      <tr>
        <td id="L2003" class="blob-num js-line-number" data-line-number="2003"></td>
        <td id="LC2003" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2004" class="blob-num js-line-number" data-line-number="2004"></td>
        <td id="LC2004" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2005" class="blob-num js-line-number" data-line-number="2005"></td>
        <td id="LC2005" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2006" class="blob-num js-line-number" data-line-number="2006"></td>
        <td id="LC2006" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2007" class="blob-num js-line-number" data-line-number="2007"></td>
        <td id="LC2007" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2008" class="blob-num js-line-number" data-line-number="2008"></td>
        <td id="LC2008" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2009" class="blob-num js-line-number" data-line-number="2009"></td>
        <td id="LC2009" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">mapDataFromRave</span>(PolarScan_t* scan, SCANMETA* meta, <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* values, <span class="pl-k">char</span>* paramStr) {</td>
      </tr>
      <tr>
        <td id="L2010" class="blob-num js-line-number" data-line-number="2010"></td>
        <td id="LC2010" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2011" class="blob-num js-line-number" data-line-number="2011"></td>
        <td id="LC2011" class="blob-code js-file-line">    PolarScanParam_t* param = <span class="pl-c1">PolarScan_getParameter</span>(scan,paramStr);</td>
      </tr>
      <tr>
        <td id="L2012" class="blob-num js-line-number" data-line-number="2012"></td>
        <td id="LC2012" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2013" class="blob-num js-line-number" data-line-number="2013"></td>
        <td id="LC2013" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L2014" class="blob-num js-line-number" data-line-number="2014"></td>
        <td id="LC2014" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzim;</td>
      </tr>
      <tr>
        <td id="L2015" class="blob-num js-line-number" data-line-number="2015"></td>
        <td id="LC2015" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2016" class="blob-num js-line-number" data-line-number="2016"></td>
        <td id="LC2016" class="blob-code js-file-line">    <span class="pl-k">int</span> nRang = <span class="pl-c1">PolarScan_getNbins</span>(scan);</td>
      </tr>
      <tr>
        <td id="L2017" class="blob-num js-line-number" data-line-number="2017"></td>
        <td id="LC2017" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzim = <span class="pl-c1">PolarScan_getNrays</span>(scan);</td>
      </tr>
      <tr>
        <td id="L2018" class="blob-num js-line-number" data-line-number="2018"></td>
        <td id="LC2018" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L2019" class="blob-num js-line-number" data-line-number="2019"></td>
        <td id="LC2019" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2020" class="blob-num js-line-number" data-line-number="2020"></td>
        <td id="LC2020" class="blob-code js-file-line">    iGlobal = <span class="pl-c1">0</span>;    </td>
      </tr>
      <tr>
        <td id="L2021" class="blob-num js-line-number" data-line-number="2021"></td>
        <td id="LC2021" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2022" class="blob-num js-line-number" data-line-number="2022"></td>
        <td id="LC2022" class="blob-code js-file-line">    <span class="pl-k">if</span> (param != <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L2023" class="blob-num js-line-number" data-line-number="2023"></td>
        <td id="LC2023" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2024" class="blob-num js-line-number" data-line-number="2024"></td>
        <td id="LC2024" class="blob-code js-file-line">        <span class="pl-c">// scan includes requested information, apparently</span></td>
      </tr>
      <tr>
        <td id="L2025" class="blob-num js-line-number" data-line-number="2025"></td>
        <td id="LC2025" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2026" class="blob-num js-line-number" data-line-number="2026"></td>
        <td id="LC2026" class="blob-code js-file-line">        meta-&gt;heig = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScan_getHeight</span>(scan);</td>
      </tr>
      <tr>
        <td id="L2027" class="blob-num js-line-number" data-line-number="2027"></td>
        <td id="LC2027" class="blob-code js-file-line">        meta-&gt;elev = (<span class="pl-k">float</span>) (<span class="pl-c1">360</span> * <span class="pl-c1">PolarScan_getElangle</span>(scan) / <span class="pl-c1">2</span> / PI );</td>
      </tr>
      <tr>
        <td id="L2028" class="blob-num js-line-number" data-line-number="2028"></td>
        <td id="LC2028" class="blob-code js-file-line">        meta-&gt;nRang = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNbins</span>(scan);</td>
      </tr>
      <tr>
        <td id="L2029" class="blob-num js-line-number" data-line-number="2029"></td>
        <td id="LC2029" class="blob-code js-file-line">        meta-&gt;nAzim = (<span class="pl-k">int</span>) <span class="pl-c1">PolarScan_getNrays</span>(scan);</td>
      </tr>
      <tr>
        <td id="L2030" class="blob-num js-line-number" data-line-number="2030"></td>
        <td id="LC2030" class="blob-code js-file-line">        meta-&gt;rangeScale = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScan_getRscale</span>(scan);</td>
      </tr>
      <tr>
        <td id="L2031" class="blob-num js-line-number" data-line-number="2031"></td>
        <td id="LC2031" class="blob-code js-file-line">        meta-&gt;azimScale = <span class="pl-c1">360</span>.<span class="pl-c1">0f</span>/meta-&gt;nAzim;   <span class="pl-c">// FIXME not sure if this always works</span></td>
      </tr>
      <tr>
        <td id="L2032" class="blob-num js-line-number" data-line-number="2032"></td>
        <td id="LC2032" class="blob-code js-file-line">        meta-&gt;valueOffset = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScanParam_getOffset</span>(param);</td>
      </tr>
      <tr>
        <td id="L2033" class="blob-num js-line-number" data-line-number="2033"></td>
        <td id="LC2033" class="blob-code js-file-line">        meta-&gt;valueScale = (<span class="pl-k">float</span>) <span class="pl-c1">PolarScanParam_getGain</span>(param);</td>
      </tr>
      <tr>
        <td id="L2034" class="blob-num js-line-number" data-line-number="2034"></td>
        <td id="LC2034" class="blob-code js-file-line">        meta-&gt;missing = (<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) <span class="pl-c1">PolarScanParam_getNodata</span>(param);</td>
      </tr>
      <tr>
        <td id="L2035" class="blob-num js-line-number" data-line-number="2035"></td>
        <td id="LC2035" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2036" class="blob-num js-line-number" data-line-number="2036"></td>
        <td id="LC2036" class="blob-code js-file-line">        <span class="pl-k">double</span> value;</td>
      </tr>
      <tr>
        <td id="L2037" class="blob-num js-line-number" data-line-number="2037"></td>
        <td id="LC2037" class="blob-code js-file-line">        <span class="pl-k">double</span>* valuePtr = &amp;value;</td>
      </tr>
      <tr>
        <td id="L2038" class="blob-num js-line-number" data-line-number="2038"></td>
        <td id="LC2038" class="blob-code js-file-line">        <span class="pl-k">unsigned</span> <span class="pl-k">char</span> valueUChar;</td>
      </tr>
      <tr>
        <td id="L2039" class="blob-num js-line-number" data-line-number="2039"></td>
        <td id="LC2039" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2040" class="blob-num js-line-number" data-line-number="2040"></td>
        <td id="LC2040" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2041" class="blob-num js-line-number" data-line-number="2041"></td>
        <td id="LC2041" class="blob-code js-file-line">        <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) {</td>
      </tr>
      <tr>
        <td id="L2042" class="blob-num js-line-number" data-line-number="2042"></td>
        <td id="LC2042" class="blob-code js-file-line">            <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L2043" class="blob-num js-line-number" data-line-number="2043"></td>
        <td id="LC2043" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2044" class="blob-num js-line-number" data-line-number="2044"></td>
        <td id="LC2044" class="blob-code js-file-line">                RaveValueType t = <span class="pl-c1">PolarScanParam_getValue</span>(param, iRang, iAzim, valuePtr);</td>
      </tr>
      <tr>
        <td id="L2045" class="blob-num js-line-number" data-line-number="2045"></td>
        <td id="LC2045" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2046" class="blob-num js-line-number" data-line-number="2046"></td>
        <td id="LC2046" class="blob-code js-file-line">                <span class="pl-k">if</span> (<span class="pl-c1">0</span> &lt;= value &amp;&amp; value &lt;=<span class="pl-c1">255</span>) {</td>
      </tr>
      <tr>
        <td id="L2047" class="blob-num js-line-number" data-line-number="2047"></td>
        <td id="LC2047" class="blob-code js-file-line">                    valueUChar = (<span class="pl-k">unsigned</span> <span class="pl-k">char</span>) value;</td>
      </tr>
      <tr>
        <td id="L2048" class="blob-num js-line-number" data-line-number="2048"></td>
        <td id="LC2048" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L2049" class="blob-num js-line-number" data-line-number="2049"></td>
        <td id="LC2049" class="blob-code js-file-line">                <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2050" class="blob-num js-line-number" data-line-number="2050"></td>
        <td id="LC2050" class="blob-code js-file-line">                    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Out of range value read at iRang = <span class="pl-c1">%d</span>, iAzim = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iRang,iAzim);</td>
      </tr>
      <tr>
        <td id="L2051" class="blob-num js-line-number" data-line-number="2051"></td>
        <td id="LC2051" class="blob-code js-file-line">                    <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2052" class="blob-num js-line-number" data-line-number="2052"></td>
        <td id="LC2052" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L2053" class="blob-num js-line-number" data-line-number="2053"></td>
        <td id="LC2053" class="blob-code js-file-line">                values[iGlobal] = valueUChar;</td>
      </tr>
      <tr>
        <td id="L2054" class="blob-num js-line-number" data-line-number="2054"></td>
        <td id="LC2054" class="blob-code js-file-line">                iGlobal++;</td>
      </tr>
      <tr>
        <td id="L2055" class="blob-num js-line-number" data-line-number="2055"></td>
        <td id="LC2055" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L2056" class="blob-num js-line-number" data-line-number="2056"></td>
        <td id="LC2056" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2057" class="blob-num js-line-number" data-line-number="2057"></td>
        <td id="LC2057" class="blob-code js-file-line">    } </td>
      </tr>
      <tr>
        <td id="L2058" class="blob-num js-line-number" data-line-number="2058"></td>
        <td id="LC2058" class="blob-code js-file-line">    <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2059" class="blob-num js-line-number" data-line-number="2059"></td>
        <td id="LC2059" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2060" class="blob-num js-line-number" data-line-number="2060"></td>
        <td id="LC2060" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Requested parameter not available<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2061" class="blob-num js-line-number" data-line-number="2061"></td>
        <td id="LC2061" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2062" class="blob-num js-line-number" data-line-number="2062"></td>
        <td id="LC2062" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2063" class="blob-num js-line-number" data-line-number="2063"></td>
        <td id="LC2063" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2064" class="blob-num js-line-number" data-line-number="2064"></td>
        <td id="LC2064" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2065" class="blob-num js-line-number" data-line-number="2065"></td>
        <td id="LC2065" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2066" class="blob-num js-line-number" data-line-number="2066"></td>
        <td id="LC2066" class="blob-code js-file-line">    <span class="pl-k">return</span> <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2067" class="blob-num js-line-number" data-line-number="2067"></td>
        <td id="LC2067" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2068" class="blob-num js-line-number" data-line-number="2068"></td>
        <td id="LC2068" class="blob-code js-file-line">} <span class="pl-c">// mapDataFromRave</span></td>
      </tr>
      <tr>
        <td id="L2069" class="blob-num js-line-number" data-line-number="2069"></td>
        <td id="LC2069" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2070" class="blob-num js-line-number" data-line-number="2070"></td>
        <td id="LC2070" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2071" class="blob-num js-line-number" data-line-number="2071"></td>
        <td id="LC2071" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2072" class="blob-num js-line-number" data-line-number="2072"></td>
        <td id="LC2072" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2073" class="blob-num js-line-number" data-line-number="2073"></td>
        <td id="LC2073" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">mapDataToRave</span>(<span class="pl-k">void</span>) {</td>
      </tr>
      <tr>
        <td id="L2074" class="blob-num js-line-number" data-line-number="2074"></td>
        <td id="LC2074" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2075" class="blob-num js-line-number" data-line-number="2075"></td>
        <td id="LC2075" class="blob-code js-file-line">    <span class="pl-c">// initialize the profile object</span></td>
      </tr>
      <tr>
        <td id="L2076" class="blob-num js-line-number" data-line-number="2076"></td>
        <td id="LC2076" class="blob-code js-file-line">    VerticalProfile_t* profileRave = <span class="pl-c1">NULL</span>;</td>
      </tr>
      <tr>
        <td id="L2077" class="blob-num js-line-number" data-line-number="2077"></td>
        <td id="LC2077" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2078" class="blob-num js-line-number" data-line-number="2078"></td>
        <td id="LC2078" class="blob-code js-file-line">    <span class="pl-c1">VerticalProfile_setLevels</span>(profileRave,nLayers);</td>
      </tr>
      <tr>
        <td id="L2079" class="blob-num js-line-number" data-line-number="2079"></td>
        <td id="LC2079" class="blob-code js-file-line">    <span class="pl-c1">VerticalProfile_setInterval</span>(profileRave,layerThickness);</td>
      </tr>
      <tr>
        <td id="L2080" class="blob-num js-line-number" data-line-number="2080"></td>
        <td id="LC2080" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2081" class="blob-num js-line-number" data-line-number="2081"></td>
        <td id="LC2081" class="blob-code js-file-line">    <span class="pl-c1">RAVE_OBJECT_RELEASE</span>(profileRave);</td>
      </tr>
      <tr>
        <td id="L2082" class="blob-num js-line-number" data-line-number="2082"></td>
        <td id="LC2082" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2083" class="blob-num js-line-number" data-line-number="2083"></td>
        <td id="LC2083" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2084" class="blob-num js-line-number" data-line-number="2084"></td>
        <td id="LC2084" class="blob-code js-file-line">    <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2085" class="blob-num js-line-number" data-line-number="2085"></td>
        <td id="LC2085" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2086" class="blob-num js-line-number" data-line-number="2086"></td>
        <td id="LC2086" class="blob-code js-file-line">}</td>
      </tr>
      <tr>
        <td id="L2087" class="blob-num js-line-number" data-line-number="2087"></td>
        <td id="LC2087" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2088" class="blob-num js-line-number" data-line-number="2088"></td>
        <td id="LC2088" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2089" class="blob-num js-line-number" data-line-number="2089"></td>
        <td id="LC2089" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2090" class="blob-num js-line-number" data-line-number="2090"></td>
        <td id="LC2090" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2091" class="blob-num js-line-number" data-line-number="2091"></td>
        <td id="LC2091" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">printGateCode</span>(<span class="pl-k">char</span>* flags, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">int</span> gateCode) {</td>
      </tr>
      <tr>
        <td id="L2092" class="blob-num js-line-number" data-line-number="2092"></td>
        <td id="LC2092" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2093" class="blob-num js-line-number" data-line-number="2093"></td>
        <td id="LC2093" class="blob-code js-file-line">    <span class="pl-c">// --------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2094" class="blob-num js-line-number" data-line-number="2094"></td>
        <td id="LC2094" class="blob-code js-file-line">    <span class="pl-c">// this function prints the integer value gateCode as  //</span></td>
      </tr>
      <tr>
        <td id="L2095" class="blob-num js-line-number" data-line-number="2095"></td>
        <td id="LC2095" class="blob-code js-file-line">    <span class="pl-c">// the equivalent sequence of bits                     //</span></td>
      </tr>
      <tr>
        <td id="L2096" class="blob-num js-line-number" data-line-number="2096"></td>
        <td id="LC2096" class="blob-code js-file-line">    <span class="pl-c">// --------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2097" class="blob-num js-line-number" data-line-number="2097"></td>
        <td id="LC2097" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2098" class="blob-num js-line-number" data-line-number="2098"></td>
        <td id="LC2098" class="blob-code js-file-line">    <span class="pl-k">int</span> iFlag;</td>
      </tr>
      <tr>
        <td id="L2099" class="blob-num js-line-number" data-line-number="2099"></td>
        <td id="LC2099" class="blob-code js-file-line">    <span class="pl-k">int</span> nFlagsNeeded;</td>
      </tr>
      <tr>
        <td id="L2100" class="blob-num js-line-number" data-line-number="2100"></td>
        <td id="LC2100" class="blob-code js-file-line">    <span class="pl-k">int</span> nFlags;</td>
      </tr>
      <tr>
        <td id="L2101" class="blob-num js-line-number" data-line-number="2101"></td>
        <td id="LC2101" class="blob-code js-file-line">    <span class="pl-k">int</span> nFlagsMax;</td>
      </tr>
      <tr>
        <td id="L2102" class="blob-num js-line-number" data-line-number="2102"></td>
        <td id="LC2102" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2103" class="blob-num js-line-number" data-line-number="2103"></td>
        <td id="LC2103" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2104" class="blob-num js-line-number" data-line-number="2104"></td>
        <td id="LC2104" class="blob-code js-file-line">    <span class="pl-k">if</span> (gateCode &lt;= <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L2105" class="blob-num js-line-number" data-line-number="2105"></td>
        <td id="LC2105" class="blob-code js-file-line">        nFlagsNeeded = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2106" class="blob-num js-line-number" data-line-number="2106"></td>
        <td id="LC2106" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2107" class="blob-num js-line-number" data-line-number="2107"></td>
        <td id="LC2107" class="blob-code js-file-line">    <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2108" class="blob-num js-line-number" data-line-number="2108"></td>
        <td id="LC2108" class="blob-code js-file-line">        nFlagsNeeded = (<span class="pl-k">int</span>) <span class="pl-c1">ceil</span>(<span class="pl-c1">log</span>(gateCode + <span class="pl-c1">1</span>)/<span class="pl-c1">log</span>(<span class="pl-c1">2</span>));</td>
      </tr>
      <tr>
        <td id="L2109" class="blob-num js-line-number" data-line-number="2109"></td>
        <td id="LC2109" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2110" class="blob-num js-line-number" data-line-number="2110"></td>
        <td id="LC2110" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2111" class="blob-num js-line-number" data-line-number="2111"></td>
        <td id="LC2111" class="blob-code js-file-line">    nFlagsMax = <span class="pl-c1">9</span>;</td>
      </tr>
      <tr>
        <td id="L2112" class="blob-num js-line-number" data-line-number="2112"></td>
        <td id="LC2112" class="blob-code js-file-line">    <span class="pl-k">if</span> (nFlagsNeeded &gt; nFlagsMax) {</td>
      </tr>
      <tr>
        <td id="L2113" class="blob-num js-line-number" data-line-number="2113"></td>
        <td id="LC2113" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>There&#39;s only space for <span class="pl-c1">%d</span> flags<span class="pl-cce">\n</span>. Aborting<span class="pl-pds">&quot;</span></span>,nFlagsMax);</td>
      </tr>
      <tr>
        <td id="L2114" class="blob-num js-line-number" data-line-number="2114"></td>
        <td id="LC2114" class="blob-code js-file-line">        <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2115" class="blob-num js-line-number" data-line-number="2115"></td>
        <td id="LC2115" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2116" class="blob-num js-line-number" data-line-number="2116"></td>
        <td id="LC2116" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2117" class="blob-num js-line-number" data-line-number="2117"></td>
        <td id="LC2117" class="blob-code js-file-line">    nFlags = nFlagsMax;</td>
      </tr>
      <tr>
        <td id="L2118" class="blob-num js-line-number" data-line-number="2118"></td>
        <td id="LC2118" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2119" class="blob-num js-line-number" data-line-number="2119"></td>
        <td id="LC2119" class="blob-code js-file-line">    <span class="pl-k">for</span> (iFlag = nFlags-<span class="pl-c1">1</span>; iFlag &gt;= <span class="pl-c1">0</span>; iFlag--) {</td>
      </tr>
      <tr>
        <td id="L2120" class="blob-num js-line-number" data-line-number="2120"></td>
        <td id="LC2120" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2121" class="blob-num js-line-number" data-line-number="2121"></td>
        <td id="LC2121" class="blob-code js-file-line">        <span class="pl-k">int</span> iFlagIsEnabled = (gateCode &amp; (<span class="pl-c1">1</span> &lt;&lt; iFlag)) &gt;&gt; iFlag;</td>
      </tr>
      <tr>
        <td id="L2122" class="blob-num js-line-number" data-line-number="2122"></td>
        <td id="LC2122" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2123" class="blob-num js-line-number" data-line-number="2123"></td>
        <td id="LC2123" class="blob-code js-file-line">        <span class="pl-k">if</span> (iFlagIsEnabled == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L2124" class="blob-num js-line-number" data-line-number="2124"></td>
        <td id="LC2124" class="blob-code js-file-line">            flags[nFlags-iFlag-<span class="pl-c1">1</span>] = <span class="pl-s"><span class="pl-pds">&#39;</span>1<span class="pl-pds">&#39;</span></span>;</td>
      </tr>
      <tr>
        <td id="L2125" class="blob-num js-line-number" data-line-number="2125"></td>
        <td id="LC2125" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2126" class="blob-num js-line-number" data-line-number="2126"></td>
        <td id="LC2126" class="blob-code js-file-line">        <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2127" class="blob-num js-line-number" data-line-number="2127"></td>
        <td id="LC2127" class="blob-code js-file-line">            flags[nFlags-iFlag-<span class="pl-c1">1</span>] = <span class="pl-s"><span class="pl-pds">&#39;</span>0<span class="pl-pds">&#39;</span></span>;</td>
      </tr>
      <tr>
        <td id="L2128" class="blob-num js-line-number" data-line-number="2128"></td>
        <td id="LC2128" class="blob-code js-file-line">        };</td>
      </tr>
      <tr>
        <td id="L2129" class="blob-num js-line-number" data-line-number="2129"></td>
        <td id="LC2129" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2130" class="blob-num js-line-number" data-line-number="2130"></td>
        <td id="LC2130" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2131" class="blob-num js-line-number" data-line-number="2131"></td>
        <td id="LC2131" class="blob-code js-file-line">    flags[nFlags] = <span class="pl-s"><span class="pl-pds">&#39;</span><span class="pl-cce">\0</span><span class="pl-pds">&#39;</span></span>;</td>
      </tr>
      <tr>
        <td id="L2132" class="blob-num js-line-number" data-line-number="2132"></td>
        <td id="LC2132" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2133" class="blob-num js-line-number" data-line-number="2133"></td>
        <td id="LC2133" class="blob-code js-file-line">    <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2134" class="blob-num js-line-number" data-line-number="2134"></td>
        <td id="LC2134" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2135" class="blob-num js-line-number" data-line-number="2135"></td>
        <td id="LC2135" class="blob-code js-file-line">} <span class="pl-c">// printGateCode</span></td>
      </tr>
      <tr>
        <td id="L2136" class="blob-num js-line-number" data-line-number="2136"></td>
        <td id="LC2136" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2137" class="blob-num js-line-number" data-line-number="2137"></td>
        <td id="LC2137" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2138" class="blob-num js-line-number" data-line-number="2138"></td>
        <td id="LC2138" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2139" class="blob-num js-line-number" data-line-number="2139"></td>
        <td id="LC2139" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2140" class="blob-num js-line-number" data-line-number="2140"></td>
        <td id="LC2140" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2141" class="blob-num js-line-number" data-line-number="2141"></td>
        <td id="LC2141" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">printMeta</span>(<span class="pl-k">const</span> SCANMETA* meta, <span class="pl-k">const</span> <span class="pl-k">char</span>* varName) {</td>
      </tr>
      <tr>
        <td id="L2142" class="blob-num js-line-number" data-line-number="2142"></td>
        <td id="LC2142" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2143" class="blob-num js-line-number" data-line-number="2143"></td>
        <td id="LC2143" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;heig = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;heig);</td>
      </tr>
      <tr>
        <td id="L2144" class="blob-num js-line-number" data-line-number="2144"></td>
        <td id="LC2144" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;elev = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;elev);</td>
      </tr>
      <tr>
        <td id="L2145" class="blob-num js-line-number" data-line-number="2145"></td>
        <td id="LC2145" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;nRang = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;nRang);</td>
      </tr>
      <tr>
        <td id="L2146" class="blob-num js-line-number" data-line-number="2146"></td>
        <td id="LC2146" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;nAzim = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;nAzim);</td>
      </tr>
      <tr>
        <td id="L2147" class="blob-num js-line-number" data-line-number="2147"></td>
        <td id="LC2147" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;rangeScale = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;rangeScale);</td>
      </tr>
      <tr>
        <td id="L2148" class="blob-num js-line-number" data-line-number="2148"></td>
        <td id="LC2148" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;azimScale = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;azimScale);</td>
      </tr>
      <tr>
        <td id="L2149" class="blob-num js-line-number" data-line-number="2149"></td>
        <td id="LC2149" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;valueOffset = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;valueOffset);</td>
      </tr>
      <tr>
        <td id="L2150" class="blob-num js-line-number" data-line-number="2150"></td>
        <td id="LC2150" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;valueScale = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;valueScale);</td>
      </tr>
      <tr>
        <td id="L2151" class="blob-num js-line-number" data-line-number="2151"></td>
        <td id="LC2151" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%s</span>-&gt;missing = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,varName,meta-&gt;missing);</td>
      </tr>
      <tr>
        <td id="L2152" class="blob-num js-line-number" data-line-number="2152"></td>
        <td id="LC2152" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2153" class="blob-num js-line-number" data-line-number="2153"></td>
        <td id="LC2153" class="blob-code js-file-line">    <span class="pl-k">return</span> <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2154" class="blob-num js-line-number" data-line-number="2154"></td>
        <td id="LC2154" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2155" class="blob-num js-line-number" data-line-number="2155"></td>
        <td id="LC2155" class="blob-code js-file-line">} <span class="pl-c">// printMeta</span></td>
      </tr>
      <tr>
        <td id="L2156" class="blob-num js-line-number" data-line-number="2156"></td>
        <td id="LC2156" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2157" class="blob-num js-line-number" data-line-number="2157"></td>
        <td id="LC2157" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2158" class="blob-num js-line-number" data-line-number="2158"></td>
        <td id="LC2158" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2159" class="blob-num js-line-number" data-line-number="2159"></td>
        <td id="LC2159" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2160" class="blob-num js-line-number" data-line-number="2160"></td>
        <td id="LC2160" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2161" class="blob-num js-line-number" data-line-number="2161"></td>
        <td id="LC2161" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">sortCells</span>(CELLPROP *cellProp, <span class="pl-k">const</span> <span class="pl-k">int</span> nCells) {</td>
      </tr>
      <tr>
        <td id="L2162" class="blob-num js-line-number" data-line-number="2162"></td>
        <td id="LC2162" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2163" class="blob-num js-line-number" data-line-number="2163"></td>
        <td id="LC2163" class="blob-code js-file-line">    <span class="pl-c">// ---------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2164" class="blob-num js-line-number" data-line-number="2164"></td>
        <td id="LC2164" class="blob-code js-file-line">    <span class="pl-c">// Sorting of the cell properties based on cell area.               //</span></td>
      </tr>
      <tr>
        <td id="L2165" class="blob-num js-line-number" data-line-number="2165"></td>
        <td id="LC2165" class="blob-code js-file-line">    <span class="pl-c">// Assume an area equal to zero for cells that are marked &#39;dropped&#39; //</span></td>
      </tr>
      <tr>
        <td id="L2166" class="blob-num js-line-number" data-line-number="2166"></td>
        <td id="LC2166" class="blob-code js-file-line">    <span class="pl-c">// ---------------------------------------------------------------- // </span></td>
      </tr>
      <tr>
        <td id="L2167" class="blob-num js-line-number" data-line-number="2167"></td>
        <td id="LC2167" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2168" class="blob-num js-line-number" data-line-number="2168"></td>
        <td id="LC2168" class="blob-code js-file-line">    <span class="pl-k">int</span> iCell;</td>
      </tr>
      <tr>
        <td id="L2169" class="blob-num js-line-number" data-line-number="2169"></td>
        <td id="LC2169" class="blob-code js-file-line">    <span class="pl-k">int</span> iCellOther;</td>
      </tr>
      <tr>
        <td id="L2170" class="blob-num js-line-number" data-line-number="2170"></td>
        <td id="LC2170" class="blob-code js-file-line">    CELLPROP tmp;</td>
      </tr>
      <tr>
        <td id="L2171" class="blob-num js-line-number" data-line-number="2171"></td>
        <td id="LC2171" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2172" class="blob-num js-line-number" data-line-number="2172"></td>
        <td id="LC2172" class="blob-code js-file-line">    <span class="pl-c">/*Sorting of data elements using straight insertion method.*/</span></td>
      </tr>
      <tr>
        <td id="L2173" class="blob-num js-line-number" data-line-number="2173"></td>
        <td id="LC2173" class="blob-code js-file-line">    <span class="pl-k">for</span> (iCell = <span class="pl-c1">1</span>; iCell &lt; nCells; iCell++) {</td>
      </tr>
      <tr>
        <td id="L2174" class="blob-num js-line-number" data-line-number="2174"></td>
        <td id="LC2174" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2175" class="blob-num js-line-number" data-line-number="2175"></td>
        <td id="LC2175" class="blob-code js-file-line">        tmp = cellProp[iCell];</td>
      </tr>
      <tr>
        <td id="L2176" class="blob-num js-line-number" data-line-number="2176"></td>
        <td id="LC2176" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2177" class="blob-num js-line-number" data-line-number="2177"></td>
        <td id="LC2177" class="blob-code js-file-line">        iCellOther = iCell - <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2178" class="blob-num js-line-number" data-line-number="2178"></td>
        <td id="LC2178" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2179" class="blob-num js-line-number" data-line-number="2179"></td>
        <td id="LC2179" class="blob-code js-file-line">        <span class="pl-k">while</span> (iCellOther &gt;= <span class="pl-c1">0</span> &amp;&amp; cellProp[iCellOther].<span class="pl-smi">nGates</span> * <span class="pl-c1">XABS</span>(cellProp[iCellOther].<span class="pl-smi">drop</span> - <span class="pl-c1">1</span>) &lt; tmp.<span class="pl-smi">nGates</span> * <span class="pl-c1">XABS</span>(tmp.<span class="pl-smi">drop</span> - <span class="pl-c1">1</span>)) {</td>
      </tr>
      <tr>
        <td id="L2180" class="blob-num js-line-number" data-line-number="2180"></td>
        <td id="LC2180" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2181" class="blob-num js-line-number" data-line-number="2181"></td>
        <td id="LC2181" class="blob-code js-file-line">            cellProp[iCellOther + <span class="pl-c1">1</span>] = cellProp[iCellOther];</td>
      </tr>
      <tr>
        <td id="L2182" class="blob-num js-line-number" data-line-number="2182"></td>
        <td id="LC2182" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2183" class="blob-num js-line-number" data-line-number="2183"></td>
        <td id="LC2183" class="blob-code js-file-line">            iCellOther--;</td>
      </tr>
      <tr>
        <td id="L2184" class="blob-num js-line-number" data-line-number="2184"></td>
        <td id="LC2184" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2185" class="blob-num js-line-number" data-line-number="2185"></td>
        <td id="LC2185" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2186" class="blob-num js-line-number" data-line-number="2186"></td>
        <td id="LC2186" class="blob-code js-file-line">        cellProp[iCellOther + <span class="pl-c1">1</span>] = tmp;</td>
      </tr>
      <tr>
        <td id="L2187" class="blob-num js-line-number" data-line-number="2187"></td>
        <td id="LC2187" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2188" class="blob-num js-line-number" data-line-number="2188"></td>
        <td id="LC2188" class="blob-code js-file-line">    } <span class="pl-c">//for iCell</span></td>
      </tr>
      <tr>
        <td id="L2189" class="blob-num js-line-number" data-line-number="2189"></td>
        <td id="LC2189" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2190" class="blob-num js-line-number" data-line-number="2190"></td>
        <td id="LC2190" class="blob-code js-file-line">    <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2191" class="blob-num js-line-number" data-line-number="2191"></td>
        <td id="LC2191" class="blob-code js-file-line">} <span class="pl-c">// sortCells</span></td>
      </tr>
      <tr>
        <td id="L2192" class="blob-num js-line-number" data-line-number="2192"></td>
        <td id="LC2192" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2193" class="blob-num js-line-number" data-line-number="2193"></td>
        <td id="LC2193" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2194" class="blob-num js-line-number" data-line-number="2194"></td>
        <td id="LC2194" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2195" class="blob-num js-line-number" data-line-number="2195"></td>
        <td id="LC2195" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2196" class="blob-num js-line-number" data-line-number="2196"></td>
        <td id="LC2196" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2197" class="blob-num js-line-number" data-line-number="2197"></td>
        <td id="LC2197" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">updateFlagFieldsInPointsArray</span>(<span class="pl-k">const</span> <span class="pl-k">float</span>* yObs, <span class="pl-k">const</span> <span class="pl-k">float</span>* yFitted, <span class="pl-k">const</span> <span class="pl-k">int</span>* includedIndex, </td>
      </tr>
      <tr>
        <td id="L2198" class="blob-num js-line-number" data-line-number="2198"></td>
        <td id="LC2198" class="blob-code js-file-line">                                   <span class="pl-k">const</span> <span class="pl-k">int</span> nPointsIncluded, <span class="pl-k">float</span>* points) {</td>
      </tr>
      <tr>
        <td id="L2199" class="blob-num js-line-number" data-line-number="2199"></td>
        <td id="LC2199" class="blob-code js-file-line">                                       </td>
      </tr>
      <tr>
        <td id="L2200" class="blob-num js-line-number" data-line-number="2200"></td>
        <td id="LC2200" class="blob-code js-file-line">    <span class="pl-c">// ----------------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2201" class="blob-num js-line-number" data-line-number="2201"></td>
        <td id="LC2201" class="blob-code js-file-line">    <span class="pl-c">// after the first svdfit to the selection of points, we want to identify gates that   //</span></td>
      </tr>
      <tr>
        <td id="L2202" class="blob-num js-line-number" data-line-number="2202"></td>
        <td id="LC2202" class="blob-code js-file-line">    <span class="pl-c">// deviate strongly from the fitted vrad value (we assume that they are outliers). So, //</span></td>
      </tr>
      <tr>
        <td id="L2203" class="blob-num js-line-number" data-line-number="2203"></td>
        <td id="LC2203" class="blob-code js-file-line">    <span class="pl-c">// this function marks the corresponding points in &#39;points&#39; by enabling the            //</span></td>
      </tr>
      <tr>
        <td id="L2204" class="blob-num js-line-number" data-line-number="2204"></td>
        <td id="LC2204" class="blob-code js-file-line">    <span class="pl-c">// &#39;flagPositionVDifMax&#39; bit in &#39;gateCode&#39;, such that outliers will be rejected during //</span></td>
      </tr>
      <tr>
        <td id="L2205" class="blob-num js-line-number" data-line-number="2205"></td>
        <td id="LC2205" class="blob-code js-file-line">    <span class="pl-c">// the second svdfit iteration                                                         //</span></td>
      </tr>
      <tr>
        <td id="L2206" class="blob-num js-line-number" data-line-number="2206"></td>
        <td id="LC2206" class="blob-code js-file-line">    <span class="pl-c">// ----------------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2207" class="blob-num js-line-number" data-line-number="2207"></td>
        <td id="LC2207" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2208" class="blob-num js-line-number" data-line-number="2208"></td>
        <td id="LC2208" class="blob-code js-file-line">    <span class="pl-k">int</span> iPointIncluded;</td>
      </tr>
      <tr>
        <td id="L2209" class="blob-num js-line-number" data-line-number="2209"></td>
        <td id="LC2209" class="blob-code js-file-line">    <span class="pl-k">int</span> iPoint;</td>
      </tr>
      <tr>
        <td id="L2210" class="blob-num js-line-number" data-line-number="2210"></td>
        <td id="LC2210" class="blob-code js-file-line">    <span class="pl-k">unsigned</span> <span class="pl-k">int</span> gateCode;</td>
      </tr>
      <tr>
        <td id="L2211" class="blob-num js-line-number" data-line-number="2211"></td>
        <td id="LC2211" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2212" class="blob-num js-line-number" data-line-number="2212"></td>
        <td id="LC2212" class="blob-code js-file-line">    <span class="pl-k">for</span> (iPointIncluded = <span class="pl-c1">0</span>; iPointIncluded &lt; nPointsIncluded; iPointIncluded++) {</td>
      </tr>
      <tr>
        <td id="L2213" class="blob-num js-line-number" data-line-number="2213"></td>
        <td id="LC2213" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2214" class="blob-num js-line-number" data-line-number="2214"></td>
        <td id="LC2214" class="blob-code js-file-line">        <span class="pl-k">float</span> absVDif = <span class="pl-c1">fabs</span>(yObs[iPointIncluded]-yFitted[iPointIncluded]);</td>
      </tr>
      <tr>
        <td id="L2215" class="blob-num js-line-number" data-line-number="2215"></td>
        <td id="LC2215" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2216" class="blob-num js-line-number" data-line-number="2216"></td>
        <td id="LC2216" class="blob-code js-file-line">        <span class="pl-k">if</span> (absVDif &gt; absVDifMax) {</td>
      </tr>
      <tr>
        <td id="L2217" class="blob-num js-line-number" data-line-number="2217"></td>
        <td id="LC2217" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2218" class="blob-num js-line-number" data-line-number="2218"></td>
        <td id="LC2218" class="blob-code js-file-line">            iPoint = includedIndex[iPointIncluded];</td>
      </tr>
      <tr>
        <td id="L2219" class="blob-num js-line-number" data-line-number="2219"></td>
        <td id="LC2219" class="blob-code js-file-line">            gateCode = (<span class="pl-k">unsigned</span> <span class="pl-k">int</span>) points[iPoint * nColsPoints + gateCodeCol];</td>
      </tr>
      <tr>
        <td id="L2220" class="blob-num js-line-number" data-line-number="2220"></td>
        <td id="LC2220" class="blob-code js-file-line">            points[iPoint * nColsPoints + gateCodeCol] = (<span class="pl-k">float</span>) (gateCode |= <span class="pl-c1">1</span>&lt;&lt;flagPositionVDifMax);</td>
      </tr>
      <tr>
        <td id="L2221" class="blob-num js-line-number" data-line-number="2221"></td>
        <td id="LC2221" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2222" class="blob-num js-line-number" data-line-number="2222"></td>
        <td id="LC2222" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2223" class="blob-num js-line-number" data-line-number="2223"></td>
        <td id="LC2223" class="blob-code js-file-line">    } </td>
      </tr>
      <tr>
        <td id="L2224" class="blob-num js-line-number" data-line-number="2224"></td>
        <td id="LC2224" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2225" class="blob-num js-line-number" data-line-number="2225"></td>
        <td id="LC2225" class="blob-code js-file-line">} <span class="pl-c">// updateFlagFieldsInPointsArray</span></td>
      </tr>
      <tr>
        <td id="L2226" class="blob-num js-line-number" data-line-number="2226"></td>
        <td id="LC2226" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2227" class="blob-num js-line-number" data-line-number="2227"></td>
        <td id="LC2227" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2228" class="blob-num js-line-number" data-line-number="2228"></td>
        <td id="LC2228" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2229" class="blob-num js-line-number" data-line-number="2229"></td>
        <td id="LC2229" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2230" class="blob-num js-line-number" data-line-number="2230"></td>
        <td id="LC2230" class="blob-code js-file-line"><span class="pl-k">static</span> <span class="pl-k">int</span> <span class="pl-en">updateMap</span>(<span class="pl-k">int</span> *cellImage, <span class="pl-k">const</span> <span class="pl-k">int</span> nGlobal, CELLPROP *cellProp, <span class="pl-k">const</span> <span class="pl-k">int</span> nCells) {</td>
      </tr>
      <tr>
        <td id="L2231" class="blob-num js-line-number" data-line-number="2231"></td>
        <td id="LC2231" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2232" class="blob-num js-line-number" data-line-number="2232"></td>
        <td id="LC2232" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2233" class="blob-num js-line-number" data-line-number="2233"></td>
        <td id="LC2233" class="blob-code js-file-line">    <span class="pl-c">// This function updates cellImage by dropping cells and reindexing the map. //</span></td>
      </tr>
      <tr>
        <td id="L2234" class="blob-num js-line-number" data-line-number="2234"></td>
        <td id="LC2234" class="blob-code js-file-line">    <span class="pl-c">// Leaving index 0 unused, will be used for assigning cell fringes           //</span></td>
      </tr>
      <tr>
        <td id="L2235" class="blob-num js-line-number" data-line-number="2235"></td>
        <td id="LC2235" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2236" class="blob-num js-line-number" data-line-number="2236"></td>
        <td id="LC2236" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2237" class="blob-num js-line-number" data-line-number="2237"></td>
        <td id="LC2237" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L2238" class="blob-num js-line-number" data-line-number="2238"></td>
        <td id="LC2238" class="blob-code js-file-line">    <span class="pl-k">int</span> iCell;</td>
      </tr>
      <tr>
        <td id="L2239" class="blob-num js-line-number" data-line-number="2239"></td>
        <td id="LC2239" class="blob-code js-file-line">    <span class="pl-k">int</span> iCellNew;</td>
      </tr>
      <tr>
        <td id="L2240" class="blob-num js-line-number" data-line-number="2240"></td>
        <td id="LC2240" class="blob-code js-file-line">    <span class="pl-k">int</span> nCellsValid;</td>
      </tr>
      <tr>
        <td id="L2241" class="blob-num js-line-number" data-line-number="2241"></td>
        <td id="LC2241" class="blob-code js-file-line">    <span class="pl-k">int</span> cellImageValue;</td>
      </tr>
      <tr>
        <td id="L2242" class="blob-num js-line-number" data-line-number="2242"></td>
        <td id="LC2242" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2243" class="blob-num js-line-number" data-line-number="2243"></td>
        <td id="LC2243" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2244" class="blob-num js-line-number" data-line-number="2244"></td>
        <td id="LC2244" class="blob-code js-file-line">    #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L2245" class="blob-num js-line-number" data-line-number="2245"></td>
        <td id="LC2245" class="blob-code js-file-line">    <span class="pl-k">int</span> minValue = cellImage[<span class="pl-c1">0</span>];</td>
      </tr>
      <tr>
        <td id="L2246" class="blob-num js-line-number" data-line-number="2246"></td>
        <td id="LC2246" class="blob-code js-file-line">    <span class="pl-k">int</span> maxValue = cellImage[<span class="pl-c1">0</span>];</td>
      </tr>
      <tr>
        <td id="L2247" class="blob-num js-line-number" data-line-number="2247"></td>
        <td id="LC2247" class="blob-code js-file-line">    <span class="pl-k">for</span> (iGlobal = <span class="pl-c1">1</span>;iGlobal &lt; nGlobal;iGlobal++) {</td>
      </tr>
      <tr>
        <td id="L2248" class="blob-num js-line-number" data-line-number="2248"></td>
        <td id="LC2248" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellImage[iGlobal] &lt; minValue) {</td>
      </tr>
      <tr>
        <td id="L2249" class="blob-num js-line-number" data-line-number="2249"></td>
        <td id="LC2249" class="blob-code js-file-line">            minValue = cellImage[iGlobal];</td>
      </tr>
      <tr>
        <td id="L2250" class="blob-num js-line-number" data-line-number="2250"></td>
        <td id="LC2250" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2251" class="blob-num js-line-number" data-line-number="2251"></td>
        <td id="LC2251" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellImage[iGlobal] &gt; maxValue) {</td>
      </tr>
      <tr>
        <td id="L2252" class="blob-num js-line-number" data-line-number="2252"></td>
        <td id="LC2252" class="blob-code js-file-line">            maxValue = cellImage[iGlobal];</td>
      </tr>
      <tr>
        <td id="L2253" class="blob-num js-line-number" data-line-number="2253"></td>
        <td id="LC2253" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2254" class="blob-num js-line-number" data-line-number="2254"></td>
        <td id="LC2254" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2255" class="blob-num js-line-number" data-line-number="2255"></td>
        <td id="LC2255" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>minimum value in cellImage array = <span class="pl-c1">%d</span>.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>, minValue);</td>
      </tr>
      <tr>
        <td id="L2256" class="blob-num js-line-number" data-line-number="2256"></td>
        <td id="LC2256" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>maximum value in cellImage array = <span class="pl-c1">%d</span>.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>, maxValue);</td>
      </tr>
      <tr>
        <td id="L2257" class="blob-num js-line-number" data-line-number="2257"></td>
        <td id="LC2257" class="blob-code js-file-line">    #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L2258" class="blob-num js-line-number" data-line-number="2258"></td>
        <td id="LC2258" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2259" class="blob-num js-line-number" data-line-number="2259"></td>
        <td id="LC2259" class="blob-code js-file-line">    nCellsValid = nCells;</td>
      </tr>
      <tr>
        <td id="L2260" class="blob-num js-line-number" data-line-number="2260"></td>
        <td id="LC2260" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2261" class="blob-num js-line-number" data-line-number="2261"></td>
        <td id="LC2261" class="blob-code js-file-line">    <span class="pl-k">for</span> (iGlobal = <span class="pl-c1">0</span>; iGlobal &lt; nGlobal; iGlobal++) {</td>
      </tr>
      <tr>
        <td id="L2262" class="blob-num js-line-number" data-line-number="2262"></td>
        <td id="LC2262" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2263" class="blob-num js-line-number" data-line-number="2263"></td>
        <td id="LC2263" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellImage[iGlobal] == -<span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L2264" class="blob-num js-line-number" data-line-number="2264"></td>
        <td id="LC2264" class="blob-code js-file-line">            <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L2265" class="blob-num js-line-number" data-line-number="2265"></td>
        <td id="LC2265" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2266" class="blob-num js-line-number" data-line-number="2266"></td>
        <td id="LC2266" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2267" class="blob-num js-line-number" data-line-number="2267"></td>
        <td id="LC2267" class="blob-code js-file-line">        cellImageValue = cellImage[iGlobal];</td>
      </tr>
      <tr>
        <td id="L2268" class="blob-num js-line-number" data-line-number="2268"></td>
        <td id="LC2268" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2269" class="blob-num js-line-number" data-line-number="2269"></td>
        <td id="LC2269" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellImageValue &gt; nCells - <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L2270" class="blob-num js-line-number" data-line-number="2270"></td>
        <td id="LC2270" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>You just asked for the properties of cell <span class="pl-c1">%d</span>, which does not exist.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>, cellImageValue);</td>
      </tr>
      <tr>
        <td id="L2271" class="blob-num js-line-number" data-line-number="2271"></td>
        <td id="LC2271" class="blob-code js-file-line">            <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L2272" class="blob-num js-line-number" data-line-number="2272"></td>
        <td id="LC2272" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2273" class="blob-num js-line-number" data-line-number="2273"></td>
        <td id="LC2273" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2274" class="blob-num js-line-number" data-line-number="2274"></td>
        <td id="LC2274" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellProp[cellImageValue].<span class="pl-smi">drop</span> == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L2275" class="blob-num js-line-number" data-line-number="2275"></td>
        <td id="LC2275" class="blob-code js-file-line">            cellImage[iGlobal] = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2276" class="blob-num js-line-number" data-line-number="2276"></td>
        <td id="LC2276" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2277" class="blob-num js-line-number" data-line-number="2277"></td>
        <td id="LC2277" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2278" class="blob-num js-line-number" data-line-number="2278"></td>
        <td id="LC2278" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2279" class="blob-num js-line-number" data-line-number="2279"></td>
        <td id="LC2279" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2280" class="blob-num js-line-number" data-line-number="2280"></td>
        <td id="LC2280" class="blob-code js-file-line">    <span class="pl-c">// sort the cells by area and determine number of valid cells</span></td>
      </tr>
      <tr>
        <td id="L2281" class="blob-num js-line-number" data-line-number="2281"></td>
        <td id="LC2281" class="blob-code js-file-line">    <span class="pl-c1">sortCells</span>(cellProp, nCells);</td>
      </tr>
      <tr>
        <td id="L2282" class="blob-num js-line-number" data-line-number="2282"></td>
        <td id="LC2282" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2283" class="blob-num js-line-number" data-line-number="2283"></td>
        <td id="LC2283" class="blob-code js-file-line">    <span class="pl-k">while</span> (nCellsValid &gt; <span class="pl-c1">0</span> &amp;&amp; cellProp[nCellsValid - <span class="pl-c1">1</span>].<span class="pl-smi">nGates</span> &lt; nGatesCellMin) {</td>
      </tr>
      <tr>
        <td id="L2284" class="blob-num js-line-number" data-line-number="2284"></td>
        <td id="LC2284" class="blob-code js-file-line">        nCellsValid--;</td>
      </tr>
      <tr>
        <td id="L2285" class="blob-num js-line-number" data-line-number="2285"></td>
        <td id="LC2285" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2286" class="blob-num js-line-number" data-line-number="2286"></td>
        <td id="LC2286" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2287" class="blob-num js-line-number" data-line-number="2287"></td>
        <td id="LC2287" class="blob-code js-file-line">    #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L2288" class="blob-num js-line-number" data-line-number="2288"></td>
        <td id="LC2288" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>nCellsValid = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,nCellsValid);</td>
      </tr>
      <tr>
        <td id="L2289" class="blob-num js-line-number" data-line-number="2289"></td>
        <td id="LC2289" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2290" class="blob-num js-line-number" data-line-number="2290"></td>
        <td id="LC2290" class="blob-code js-file-line">    #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L2291" class="blob-num js-line-number" data-line-number="2291"></td>
        <td id="LC2291" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2292" class="blob-num js-line-number" data-line-number="2292"></td>
        <td id="LC2292" class="blob-code js-file-line">    <span class="pl-c">// replace the values in cellImage with newly calculated index values:</span></td>
      </tr>
      <tr>
        <td id="L2293" class="blob-num js-line-number" data-line-number="2293"></td>
        <td id="LC2293" class="blob-code js-file-line">    <span class="pl-k">for</span> (iCell = <span class="pl-c1">0</span>; iCell &lt; nCells; iCell++) {</td>
      </tr>
      <tr>
        <td id="L2294" class="blob-num js-line-number" data-line-number="2294"></td>
        <td id="LC2294" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2295" class="blob-num js-line-number" data-line-number="2295"></td>
        <td id="LC2295" class="blob-code js-file-line">        <span class="pl-k">if</span> (iCell &lt; nCellsValid) {</td>
      </tr>
      <tr>
        <td id="L2296" class="blob-num js-line-number" data-line-number="2296"></td>
        <td id="LC2296" class="blob-code js-file-line">            iCellNew = -<span class="pl-c1">1</span> * (iCell + <span class="pl-c1">2</span> + <span class="pl-c1">100</span>);</td>
      </tr>
      <tr>
        <td id="L2297" class="blob-num js-line-number" data-line-number="2297"></td>
        <td id="LC2297" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2298" class="blob-num js-line-number" data-line-number="2298"></td>
        <td id="LC2298" class="blob-code js-file-line">        <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2299" class="blob-num js-line-number" data-line-number="2299"></td>
        <td id="LC2299" class="blob-code js-file-line">            iCellNew = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2300" class="blob-num js-line-number" data-line-number="2300"></td>
        <td id="LC2300" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2301" class="blob-num js-line-number" data-line-number="2301"></td>
        <td id="LC2301" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2302" class="blob-num js-line-number" data-line-number="2302"></td>
        <td id="LC2302" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L2303" class="blob-num js-line-number" data-line-number="2303"></td>
        <td id="LC2303" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>before: cellProp[<span class="pl-c1">%d</span>].index = <span class="pl-c1">%d</span>.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iCell,cellProp[iCell].<span class="pl-smi">index</span>);</td>
      </tr>
      <tr>
        <td id="L2304" class="blob-num js-line-number" data-line-number="2304"></td>
        <td id="LC2304" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>before: cellProp[<span class="pl-c1">%d</span>].nGates = <span class="pl-c1">%d</span>.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iCell,cellProp[iCell].<span class="pl-smi">nGates</span>);</td>
      </tr>
      <tr>
        <td id="L2305" class="blob-num js-line-number" data-line-number="2305"></td>
        <td id="LC2305" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>before: iCell = <span class="pl-c1">%d</span>.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iCell);</td>
      </tr>
      <tr>
        <td id="L2306" class="blob-num js-line-number" data-line-number="2306"></td>
        <td id="LC2306" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>before: iCellNew = <span class="pl-c1">%d</span>.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iCellNew);</td>
      </tr>
      <tr>
        <td id="L2307" class="blob-num js-line-number" data-line-number="2307"></td>
        <td id="LC2307" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2308" class="blob-num js-line-number" data-line-number="2308"></td>
        <td id="LC2308" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L2309" class="blob-num js-line-number" data-line-number="2309"></td>
        <td id="LC2309" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2310" class="blob-num js-line-number" data-line-number="2310"></td>
        <td id="LC2310" class="blob-code js-file-line">        <span class="pl-k">for</span> (iGlobal = <span class="pl-c1">0</span>; iGlobal &lt; nGlobal; iGlobal++) {</td>
      </tr>
      <tr>
        <td id="L2311" class="blob-num js-line-number" data-line-number="2311"></td>
        <td id="LC2311" class="blob-code js-file-line">            <span class="pl-k">if</span> (cellImage[iGlobal] == cellProp[iCell].<span class="pl-smi">index</span>) {</td>
      </tr>
      <tr>
        <td id="L2312" class="blob-num js-line-number" data-line-number="2312"></td>
        <td id="LC2312" class="blob-code js-file-line">                cellImage[iGlobal] = iCellNew;</td>
      </tr>
      <tr>
        <td id="L2313" class="blob-num js-line-number" data-line-number="2313"></td>
        <td id="LC2313" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L2314" class="blob-num js-line-number" data-line-number="2314"></td>
        <td id="LC2314" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2315" class="blob-num js-line-number" data-line-number="2315"></td>
        <td id="LC2315" class="blob-code js-file-line">        <span class="pl-c">// have the indices in cellProp match the re-numbering</span></td>
      </tr>
      <tr>
        <td id="L2316" class="blob-num js-line-number" data-line-number="2316"></td>
        <td id="LC2316" class="blob-code js-file-line">        cellProp[iCell].<span class="pl-smi">index</span> = iCellNew;</td>
      </tr>
      <tr>
        <td id="L2317" class="blob-num js-line-number" data-line-number="2317"></td>
        <td id="LC2317" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2318" class="blob-num js-line-number" data-line-number="2318"></td>
        <td id="LC2318" class="blob-code js-file-line">    } <span class="pl-c">// (iCell = 0; iCell &lt; nCells; iCell++)</span></td>
      </tr>
      <tr>
        <td id="L2319" class="blob-num js-line-number" data-line-number="2319"></td>
        <td id="LC2319" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2320" class="blob-num js-line-number" data-line-number="2320"></td>
        <td id="LC2320" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2321" class="blob-num js-line-number" data-line-number="2321"></td>
        <td id="LC2321" class="blob-code js-file-line">    <span class="pl-c">// once you&#39;ve re-numbered everything, flip the sign back and</span></td>
      </tr>
      <tr>
        <td id="L2322" class="blob-num js-line-number" data-line-number="2322"></td>
        <td id="LC2322" class="blob-code js-file-line">    <span class="pl-c">// remove the offset of 100...</span></td>
      </tr>
      <tr>
        <td id="L2323" class="blob-num js-line-number" data-line-number="2323"></td>
        <td id="LC2323" class="blob-code js-file-line">    <span class="pl-k">for</span> (iGlobal = <span class="pl-c1">0</span>; iGlobal &lt; nGlobal; iGlobal++) {</td>
      </tr>
      <tr>
        <td id="L2324" class="blob-num js-line-number" data-line-number="2324"></td>
        <td id="LC2324" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellImage[iGlobal] == -<span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L2325" class="blob-num js-line-number" data-line-number="2325"></td>
        <td id="LC2325" class="blob-code js-file-line">            <span class="pl-c">// do nothing</span></td>
      </tr>
      <tr>
        <td id="L2326" class="blob-num js-line-number" data-line-number="2326"></td>
        <td id="LC2326" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2327" class="blob-num js-line-number" data-line-number="2327"></td>
        <td id="LC2327" class="blob-code js-file-line">        <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2328" class="blob-num js-line-number" data-line-number="2328"></td>
        <td id="LC2328" class="blob-code js-file-line">            cellImage[iGlobal] = (-<span class="pl-c1">1</span> * cellImage[iGlobal]) - <span class="pl-c1">100</span>;</td>
      </tr>
      <tr>
        <td id="L2329" class="blob-num js-line-number" data-line-number="2329"></td>
        <td id="LC2329" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2330" class="blob-num js-line-number" data-line-number="2330"></td>
        <td id="LC2330" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2331" class="blob-num js-line-number" data-line-number="2331"></td>
        <td id="LC2331" class="blob-code js-file-line">    <span class="pl-c">// ...and make sure the indices in cellProp match that change</span></td>
      </tr>
      <tr>
        <td id="L2332" class="blob-num js-line-number" data-line-number="2332"></td>
        <td id="LC2332" class="blob-code js-file-line">    <span class="pl-k">for</span> (iCell = <span class="pl-c1">0</span>; iCell &lt; nCells; iCell++) {</td>
      </tr>
      <tr>
        <td id="L2333" class="blob-num js-line-number" data-line-number="2333"></td>
        <td id="LC2333" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2334" class="blob-num js-line-number" data-line-number="2334"></td>
        <td id="LC2334" class="blob-code js-file-line">        <span class="pl-k">if</span> (cellProp[iCell].<span class="pl-smi">index</span> == -<span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L2335" class="blob-num js-line-number" data-line-number="2335"></td>
        <td id="LC2335" class="blob-code js-file-line">            <span class="pl-c">// do nothing</span></td>
      </tr>
      <tr>
        <td id="L2336" class="blob-num js-line-number" data-line-number="2336"></td>
        <td id="LC2336" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2337" class="blob-num js-line-number" data-line-number="2337"></td>
        <td id="LC2337" class="blob-code js-file-line">        <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2338" class="blob-num js-line-number" data-line-number="2338"></td>
        <td id="LC2338" class="blob-code js-file-line">            cellProp[iCell].<span class="pl-smi">index</span> = (-<span class="pl-c1">1</span> * cellProp[iCell].<span class="pl-smi">index</span>) - <span class="pl-c1">100</span>;</td>
      </tr>
      <tr>
        <td id="L2339" class="blob-num js-line-number" data-line-number="2339"></td>
        <td id="LC2339" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2340" class="blob-num js-line-number" data-line-number="2340"></td>
        <td id="LC2340" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2341" class="blob-num js-line-number" data-line-number="2341"></td>
        <td id="LC2341" class="blob-code js-file-line">        #<span class="pl-k">ifdef</span> FPRINTFON</td>
      </tr>
      <tr>
        <td id="L2342" class="blob-num js-line-number" data-line-number="2342"></td>
        <td id="LC2342" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>after: cellProp[<span class="pl-c1">%d</span>].index = <span class="pl-c1">%d</span>.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iCell,cellProp[iCell].<span class="pl-smi">index</span>);</td>
      </tr>
      <tr>
        <td id="L2343" class="blob-num js-line-number" data-line-number="2343"></td>
        <td id="LC2343" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>after: cellProp[<span class="pl-c1">%d</span>].nGates = <span class="pl-c1">%d</span>.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iCell,cellProp[iCell].<span class="pl-smi">nGates</span>);</td>
      </tr>
      <tr>
        <td id="L2344" class="blob-num js-line-number" data-line-number="2344"></td>
        <td id="LC2344" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2345" class="blob-num js-line-number" data-line-number="2345"></td>
        <td id="LC2345" class="blob-code js-file-line">        #<span class="pl-k">endif</span></td>
      </tr>
      <tr>
        <td id="L2346" class="blob-num js-line-number" data-line-number="2346"></td>
        <td id="LC2346" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2347" class="blob-num js-line-number" data-line-number="2347"></td>
        <td id="LC2347" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2348" class="blob-num js-line-number" data-line-number="2348"></td>
        <td id="LC2348" class="blob-code js-file-line">    <span class="pl-k">return</span> nCellsValid;</td>
      </tr>
      <tr>
        <td id="L2349" class="blob-num js-line-number" data-line-number="2349"></td>
        <td id="LC2349" class="blob-code js-file-line">} <span class="pl-c">// updateMap</span></td>
      </tr>
      <tr>
        <td id="L2350" class="blob-num js-line-number" data-line-number="2350"></td>
        <td id="LC2350" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2351" class="blob-num js-line-number" data-line-number="2351"></td>
        <td id="LC2351" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2352" class="blob-num js-line-number" data-line-number="2352"></td>
        <td id="LC2352" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2353" class="blob-num js-line-number" data-line-number="2353"></td>
        <td id="LC2353" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2354" class="blob-num js-line-number" data-line-number="2354"></td>
        <td id="LC2354" class="blob-code js-file-line"><span class="pl-k">void</span> <span class="pl-en">vol2birdCalcProfiles</span>() {</td>
      </tr>
      <tr>
        <td id="L2355" class="blob-num js-line-number" data-line-number="2355"></td>
        <td id="LC2355" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2356" class="blob-num js-line-number" data-line-number="2356"></td>
        <td id="LC2356" class="blob-code js-file-line">    <span class="pl-k">int</span> nPasses;</td>
      </tr>
      <tr>
        <td id="L2357" class="blob-num js-line-number" data-line-number="2357"></td>
        <td id="LC2357" class="blob-code js-file-line">    <span class="pl-k">int</span> iPoint;</td>
      </tr>
      <tr>
        <td id="L2358" class="blob-num js-line-number" data-line-number="2358"></td>
        <td id="LC2358" class="blob-code js-file-line">    <span class="pl-k">int</span> iLayer;</td>
      </tr>
      <tr>
        <td id="L2359" class="blob-num js-line-number" data-line-number="2359"></td>
        <td id="LC2359" class="blob-code js-file-line">    <span class="pl-k">int</span> iPass;</td>
      </tr>
      <tr>
        <td id="L2360" class="blob-num js-line-number" data-line-number="2360"></td>
        <td id="LC2360" class="blob-code js-file-line">    <span class="pl-k">int</span> iProfileType;</td>
      </tr>
      <tr>
        <td id="L2361" class="blob-num js-line-number" data-line-number="2361"></td>
        <td id="LC2361" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2362" class="blob-num js-line-number" data-line-number="2362"></td>
        <td id="LC2362" class="blob-code js-file-line">    <span class="pl-k">if</span> (initializationSuccessful==<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L2363" class="blob-num js-line-number" data-line-number="2363"></td>
        <td id="LC2363" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>You need to initialize vol2bird before you can use it. Aborting.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2364" class="blob-num js-line-number" data-line-number="2364"></td>
        <td id="LC2364" class="blob-code js-file-line">        <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2365" class="blob-num js-line-number" data-line-number="2365"></td>
        <td id="LC2365" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2366" class="blob-num js-line-number" data-line-number="2366"></td>
        <td id="LC2366" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2367" class="blob-num js-line-number" data-line-number="2367"></td>
        <td id="LC2367" class="blob-code js-file-line">    <span class="pl-c">// calculate the profiles in reverse order, because you need the result </span></td>
      </tr>
      <tr>
        <td id="L2368" class="blob-num js-line-number" data-line-number="2368"></td>
        <td id="LC2368" class="blob-code js-file-line">    <span class="pl-c">// of iProfileType == 3 in order to check whether chi &lt; stdDevMinBird  </span></td>
      </tr>
      <tr>
        <td id="L2369" class="blob-num js-line-number" data-line-number="2369"></td>
        <td id="LC2369" class="blob-code js-file-line">    <span class="pl-c">// when calculating iProfileType == 1</span></td>
      </tr>
      <tr>
        <td id="L2370" class="blob-num js-line-number" data-line-number="2370"></td>
        <td id="LC2370" class="blob-code js-file-line">    <span class="pl-k">for</span> (iProfileType = nProfileTypes; iProfileType &gt; <span class="pl-c1">0</span>; iProfileType--) {</td>
      </tr>
      <tr>
        <td id="L2371" class="blob-num js-line-number" data-line-number="2371"></td>
        <td id="LC2371" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2372" class="blob-num js-line-number" data-line-number="2372"></td>
        <td id="LC2372" class="blob-code js-file-line">        <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2373" class="blob-num js-line-number" data-line-number="2373"></td>
        <td id="LC2373" class="blob-code js-file-line">        <span class="pl-c">//                        prepare the profile                    //</span></td>
      </tr>
      <tr>
        <td id="L2374" class="blob-num js-line-number" data-line-number="2374"></td>
        <td id="LC2374" class="blob-code js-file-line">        <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2375" class="blob-num js-line-number" data-line-number="2375"></td>
        <td id="LC2375" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2376" class="blob-num js-line-number" data-line-number="2376"></td>
        <td id="LC2376" class="blob-code js-file-line">        iProfileTypeLast = iProfileType;</td>
      </tr>
      <tr>
        <td id="L2377" class="blob-num js-line-number" data-line-number="2377"></td>
        <td id="LC2377" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2378" class="blob-num js-line-number" data-line-number="2378"></td>
        <td id="LC2378" class="blob-code js-file-line">        <span class="pl-c">// if the user does not require fitting a model to the observed </span></td>
      </tr>
      <tr>
        <td id="L2379" class="blob-num js-line-number" data-line-number="2379"></td>
        <td id="LC2379" class="blob-code js-file-line">        <span class="pl-c">// vrad values, one pass is enough </span></td>
      </tr>
      <tr>
        <td id="L2380" class="blob-num js-line-number" data-line-number="2380"></td>
        <td id="LC2380" class="blob-code js-file-line">        <span class="pl-k">if</span> (fitVrad == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L2381" class="blob-num js-line-number" data-line-number="2381"></td>
        <td id="LC2381" class="blob-code js-file-line">            nPasses = <span class="pl-c1">2</span>;</td>
      </tr>
      <tr>
        <td id="L2382" class="blob-num js-line-number" data-line-number="2382"></td>
        <td id="LC2382" class="blob-code js-file-line">        } </td>
      </tr>
      <tr>
        <td id="L2383" class="blob-num js-line-number" data-line-number="2383"></td>
        <td id="LC2383" class="blob-code js-file-line">        <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2384" class="blob-num js-line-number" data-line-number="2384"></td>
        <td id="LC2384" class="blob-code js-file-line">            nPasses = <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2385" class="blob-num js-line-number" data-line-number="2385"></td>
        <td id="LC2385" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2386" class="blob-num js-line-number" data-line-number="2386"></td>
        <td id="LC2386" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2387" class="blob-num js-line-number" data-line-number="2387"></td>
        <td id="LC2387" class="blob-code js-file-line">        <span class="pl-c">// reset the flagPositionVDifMax bit before calculating each profile</span></td>
      </tr>
      <tr>
        <td id="L2388" class="blob-num js-line-number" data-line-number="2388"></td>
        <td id="LC2388" class="blob-code js-file-line">        <span class="pl-k">for</span> (iPoint = <span class="pl-c1">0</span>; iPoint &lt; nRowsPoints; iPoint++) {</td>
      </tr>
      <tr>
        <td id="L2389" class="blob-num js-line-number" data-line-number="2389"></td>
        <td id="LC2389" class="blob-code js-file-line">            <span class="pl-k">unsigned</span> <span class="pl-k">int</span> gateCode = (<span class="pl-k">unsigned</span> <span class="pl-k">int</span>) points[iPoint * nColsPoints + gateCodeCol];</td>
      </tr>
      <tr>
        <td id="L2390" class="blob-num js-line-number" data-line-number="2390"></td>
        <td id="LC2390" class="blob-code js-file-line">            points[iPoint * nColsPoints + gateCodeCol] = (<span class="pl-k">float</span>) (gateCode &amp;= ~(<span class="pl-c1">1</span>&lt;&lt;flagPositionVDifMax));</td>
      </tr>
      <tr>
        <td id="L2391" class="blob-num js-line-number" data-line-number="2391"></td>
        <td id="LC2391" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2392" class="blob-num js-line-number" data-line-number="2392"></td>
        <td id="LC2392" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2393" class="blob-num js-line-number" data-line-number="2393"></td>
        <td id="LC2393" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2394" class="blob-num js-line-number" data-line-number="2394"></td>
        <td id="LC2394" class="blob-code js-file-line">        <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L2395" class="blob-num js-line-number" data-line-number="2395"></td>
        <td id="LC2395" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2396" class="blob-num js-line-number" data-line-number="2396"></td>
        <td id="LC2396" class="blob-code js-file-line">            <span class="pl-c">// this variable is needed just outside of the iPass loop below </span></td>
      </tr>
      <tr>
        <td id="L2397" class="blob-num js-line-number" data-line-number="2397"></td>
        <td id="LC2397" class="blob-code js-file-line">            <span class="pl-k">float</span> chi = NAN;</td>
      </tr>
      <tr>
        <td id="L2398" class="blob-num js-line-number" data-line-number="2398"></td>
        <td id="LC2398" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2399" class="blob-num js-line-number" data-line-number="2399"></td>
        <td id="LC2399" class="blob-code js-file-line">            <span class="pl-k">for</span> (iPass = <span class="pl-c1">0</span>; iPass &lt; nPasses; iPass++) {</td>
      </tr>
      <tr>
        <td id="L2400" class="blob-num js-line-number" data-line-number="2400"></td>
        <td id="LC2400" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2401" class="blob-num js-line-number" data-line-number="2401"></td>
        <td id="LC2401" class="blob-code js-file-line">                <span class="pl-k">const</span> <span class="pl-k">int</span> iPointFrom = indexFrom[iLayer];</td>
      </tr>
      <tr>
        <td id="L2402" class="blob-num js-line-number" data-line-number="2402"></td>
        <td id="LC2402" class="blob-code js-file-line">                <span class="pl-k">const</span> <span class="pl-k">int</span> nPointsLayer = nPointsWritten[iLayer];</td>
      </tr>
      <tr>
        <td id="L2403" class="blob-num js-line-number" data-line-number="2403"></td>
        <td id="LC2403" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2404" class="blob-num js-line-number" data-line-number="2404"></td>
        <td id="LC2404" class="blob-code js-file-line">                <span class="pl-k">int</span> iPointLayer;</td>
      </tr>
      <tr>
        <td id="L2405" class="blob-num js-line-number" data-line-number="2405"></td>
        <td id="LC2405" class="blob-code js-file-line">                <span class="pl-k">int</span> iPointIncluded;</td>
      </tr>
      <tr>
        <td id="L2406" class="blob-num js-line-number" data-line-number="2406"></td>
        <td id="LC2406" class="blob-code js-file-line">                <span class="pl-k">int</span> nPointsIncluded;</td>
      </tr>
      <tr>
        <td id="L2407" class="blob-num js-line-number" data-line-number="2407"></td>
        <td id="LC2407" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2408" class="blob-num js-line-number" data-line-number="2408"></td>
        <td id="LC2408" class="blob-code js-file-line">                <span class="pl-k">float</span> parameterVector[] = {NAN,NAN,NAN};</td>
      </tr>
      <tr>
        <td id="L2409" class="blob-num js-line-number" data-line-number="2409"></td>
        <td id="LC2409" class="blob-code js-file-line">                <span class="pl-k">float</span> avar[] = {NAN,NAN,NAN};</td>
      </tr>
      <tr>
        <td id="L2410" class="blob-num js-line-number" data-line-number="2410"></td>
        <td id="LC2410" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2411" class="blob-num js-line-number" data-line-number="2411"></td>
        <td id="LC2411" class="blob-code js-file-line">                <span class="pl-k">float</span>* pointsSelection = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">float</span>) * nPointsLayer * nDims);</td>
      </tr>
      <tr>
        <td id="L2412" class="blob-num js-line-number" data-line-number="2412"></td>
        <td id="LC2412" class="blob-code js-file-line">                <span class="pl-k">float</span>* yObs = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">float</span>) * nPointsLayer);</td>
      </tr>
      <tr>
        <td id="L2413" class="blob-num js-line-number" data-line-number="2413"></td>
        <td id="LC2413" class="blob-code js-file-line">                <span class="pl-k">float</span>* yFitted = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">float</span>) * nPointsLayer);</td>
      </tr>
      <tr>
        <td id="L2414" class="blob-num js-line-number" data-line-number="2414"></td>
        <td id="LC2414" class="blob-code js-file-line">                <span class="pl-k">int</span>* includedIndex = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nPointsLayer);</td>
      </tr>
      <tr>
        <td id="L2415" class="blob-num js-line-number" data-line-number="2415"></td>
        <td id="LC2415" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2416" class="blob-num js-line-number" data-line-number="2416"></td>
        <td id="LC2416" class="blob-code js-file-line">                <span class="pl-k">float</span> dbzValue = NAN;</td>
      </tr>
      <tr>
        <td id="L2417" class="blob-num js-line-number" data-line-number="2417"></td>
        <td id="LC2417" class="blob-code js-file-line">                <span class="pl-k">float</span> undbzValue = NAN;</td>
      </tr>
      <tr>
        <td id="L2418" class="blob-num js-line-number" data-line-number="2418"></td>
        <td id="LC2418" class="blob-code js-file-line">                <span class="pl-k">double</span> undbzSum = <span class="pl-c1">0.0</span>;</td>
      </tr>
      <tr>
        <td id="L2419" class="blob-num js-line-number" data-line-number="2419"></td>
        <td id="LC2419" class="blob-code js-file-line">                <span class="pl-k">float</span> undbzAvg = NAN;</td>
      </tr>
      <tr>
        <td id="L2420" class="blob-num js-line-number" data-line-number="2420"></td>
        <td id="LC2420" class="blob-code js-file-line">                <span class="pl-k">float</span> dbzAvg = NAN;</td>
      </tr>
      <tr>
        <td id="L2421" class="blob-num js-line-number" data-line-number="2421"></td>
        <td id="LC2421" class="blob-code js-file-line">                <span class="pl-k">float</span> birdDensity = NAN;</td>
      </tr>
      <tr>
        <td id="L2422" class="blob-num js-line-number" data-line-number="2422"></td>
        <td id="LC2422" class="blob-code js-file-line">                <span class="pl-k">float</span> reflectivity = NAN;</td>
      </tr>
      <tr>
        <td id="L2423" class="blob-num js-line-number" data-line-number="2423"></td>
        <td id="LC2423" class="blob-code js-file-line">                <span class="pl-k">int</span> hasGap = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L2424" class="blob-num js-line-number" data-line-number="2424"></td>
        <td id="LC2424" class="blob-code js-file-line">                <span class="pl-k">float</span> chisq = NAN;</td>
      </tr>
      <tr>
        <td id="L2425" class="blob-num js-line-number" data-line-number="2425"></td>
        <td id="LC2425" class="blob-code js-file-line">                <span class="pl-k">float</span> hSpeed = NAN;</td>
      </tr>
      <tr>
        <td id="L2426" class="blob-num js-line-number" data-line-number="2426"></td>
        <td id="LC2426" class="blob-code js-file-line">                <span class="pl-k">float</span> hDir = NAN;</td>
      </tr>
      <tr>
        <td id="L2427" class="blob-num js-line-number" data-line-number="2427"></td>
        <td id="LC2427" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2428" class="blob-num js-line-number" data-line-number="2428"></td>
        <td id="LC2428" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2429" class="blob-num js-line-number" data-line-number="2429"></td>
        <td id="LC2429" class="blob-code js-file-line">                <span class="pl-k">for</span> (iPointLayer = <span class="pl-c1">0</span>; iPointLayer &lt; nPointsLayer; iPointLayer++) {</td>
      </tr>
      <tr>
        <td id="L2430" class="blob-num js-line-number" data-line-number="2430"></td>
        <td id="LC2430" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2431" class="blob-num js-line-number" data-line-number="2431"></td>
        <td id="LC2431" class="blob-code js-file-line">                    pointsSelection[iPointLayer * nDims + <span class="pl-c1">0</span>] = <span class="pl-c1">0</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L2432" class="blob-num js-line-number" data-line-number="2432"></td>
        <td id="LC2432" class="blob-code js-file-line">                    pointsSelection[iPointLayer * nDims + <span class="pl-c1">1</span>] = <span class="pl-c1">0</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L2433" class="blob-num js-line-number" data-line-number="2433"></td>
        <td id="LC2433" class="blob-code js-file-line">                    </td>
      </tr>
      <tr>
        <td id="L2434" class="blob-num js-line-number" data-line-number="2434"></td>
        <td id="LC2434" class="blob-code js-file-line">                    yObs[iPointLayer] = <span class="pl-c1">0</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L2435" class="blob-num js-line-number" data-line-number="2435"></td>
        <td id="LC2435" class="blob-code js-file-line">                    yFitted[iPointLayer] = <span class="pl-c1">0</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L2436" class="blob-num js-line-number" data-line-number="2436"></td>
        <td id="LC2436" class="blob-code js-file-line">                    </td>
      </tr>
      <tr>
        <td id="L2437" class="blob-num js-line-number" data-line-number="2437"></td>
        <td id="LC2437" class="blob-code js-file-line">                    includedIndex[iPointLayer] = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2438" class="blob-num js-line-number" data-line-number="2438"></td>
        <td id="LC2438" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2439" class="blob-num js-line-number" data-line-number="2439"></td>
        <td id="LC2439" class="blob-code js-file-line">                };</td>
      </tr>
      <tr>
        <td id="L2440" class="blob-num js-line-number" data-line-number="2440"></td>
        <td id="LC2440" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2441" class="blob-num js-line-number" data-line-number="2441"></td>
        <td id="LC2441" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">0</span>] = iLayer * layerThickness;</td>
      </tr>
      <tr>
        <td id="L2442" class="blob-num js-line-number" data-line-number="2442"></td>
        <td id="LC2442" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">1</span>] = (iLayer + <span class="pl-c1">1</span>) * layerThickness;</td>
      </tr>
      <tr>
        <td id="L2443" class="blob-num js-line-number" data-line-number="2443"></td>
        <td id="LC2443" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">2</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2444" class="blob-num js-line-number" data-line-number="2444"></td>
        <td id="LC2444" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">3</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2445" class="blob-num js-line-number" data-line-number="2445"></td>
        <td id="LC2445" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">4</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2446" class="blob-num js-line-number" data-line-number="2446"></td>
        <td id="LC2446" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">5</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2447" class="blob-num js-line-number" data-line-number="2447"></td>
        <td id="LC2447" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">6</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2448" class="blob-num js-line-number" data-line-number="2448"></td>
        <td id="LC2448" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">7</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2449" class="blob-num js-line-number" data-line-number="2449"></td>
        <td id="LC2449" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">8</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2450" class="blob-num js-line-number" data-line-number="2450"></td>
        <td id="LC2450" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">9</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2451" class="blob-num js-line-number" data-line-number="2451"></td>
        <td id="LC2451" class="blob-code js-file-line">                profile[iLayer*nColsProfile + <span class="pl-c1">10</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2452" class="blob-num js-line-number" data-line-number="2452"></td>
        <td id="LC2452" class="blob-code js-file-line">                profile[iLayer*nColsProfile + <span class="pl-c1">11</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2453" class="blob-num js-line-number" data-line-number="2453"></td>
        <td id="LC2453" class="blob-code js-file-line">                profile[iLayer*nColsProfile + <span class="pl-c1">12</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2454" class="blob-num js-line-number" data-line-number="2454"></td>
        <td id="LC2454" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2455" class="blob-num js-line-number" data-line-number="2455"></td>
        <td id="LC2455" class="blob-code js-file-line">                iPointIncluded = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2456" class="blob-num js-line-number" data-line-number="2456"></td>
        <td id="LC2456" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2457" class="blob-num js-line-number" data-line-number="2457"></td>
        <td id="LC2457" class="blob-code js-file-line">                <span class="pl-k">for</span> (iPointLayer = iPointFrom; iPointLayer &lt; iPointFrom + nPointsLayer; iPointLayer++) {</td>
      </tr>
      <tr>
        <td id="L2458" class="blob-num js-line-number" data-line-number="2458"></td>
        <td id="LC2458" class="blob-code js-file-line">                    </td>
      </tr>
      <tr>
        <td id="L2459" class="blob-num js-line-number" data-line-number="2459"></td>
        <td id="LC2459" class="blob-code js-file-line">                    <span class="pl-k">unsigned</span> <span class="pl-k">int</span> gateCode = (<span class="pl-k">unsigned</span> <span class="pl-k">int</span>) points[iPointLayer * nColsPoints + gateCodeCol];</td>
      </tr>
      <tr>
        <td id="L2460" class="blob-num js-line-number" data-line-number="2460"></td>
        <td id="LC2460" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2461" class="blob-num js-line-number" data-line-number="2461"></td>
        <td id="LC2461" class="blob-code js-file-line">                    <span class="pl-k">if</span> (<span class="pl-c1">includeGate</span>(iProfileType,gateCode) == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L2462" class="blob-num js-line-number" data-line-number="2462"></td>
        <td id="LC2462" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2463" class="blob-num js-line-number" data-line-number="2463"></td>
        <td id="LC2463" class="blob-code js-file-line">                        <span class="pl-c">// copy azimuth angle from the &#39;points&#39; array</span></td>
      </tr>
      <tr>
        <td id="L2464" class="blob-num js-line-number" data-line-number="2464"></td>
        <td id="LC2464" class="blob-code js-file-line">                        pointsSelection[iPointIncluded * nDims + <span class="pl-c1">0</span>] = points[iPointLayer * nColsPoints + azimAngleCol];</td>
      </tr>
      <tr>
        <td id="L2465" class="blob-num js-line-number" data-line-number="2465"></td>
        <td id="LC2465" class="blob-code js-file-line">                        <span class="pl-c">// copy elevation angle from the &#39;points&#39; array</span></td>
      </tr>
      <tr>
        <td id="L2466" class="blob-num js-line-number" data-line-number="2466"></td>
        <td id="LC2466" class="blob-code js-file-line">                        pointsSelection[iPointIncluded * nDims + <span class="pl-c1">1</span>] = points[iPointLayer * nColsPoints + elevAngleCol];</td>
      </tr>
      <tr>
        <td id="L2467" class="blob-num js-line-number" data-line-number="2467"></td>
        <td id="LC2467" class="blob-code js-file-line">                        <span class="pl-c">// get the dbz value at this [azimuth, elevation] </span></td>
      </tr>
      <tr>
        <td id="L2468" class="blob-num js-line-number" data-line-number="2468"></td>
        <td id="LC2468" class="blob-code js-file-line">                        dbzValue = points[iPointLayer * nColsPoints + dbzValueCol];</td>
      </tr>
      <tr>
        <td id="L2469" class="blob-num js-line-number" data-line-number="2469"></td>
        <td id="LC2469" class="blob-code js-file-line">                        <span class="pl-c">// convert from dB scale to linear scale </span></td>
      </tr>
      <tr>
        <td id="L2470" class="blob-num js-line-number" data-line-number="2470"></td>
        <td id="LC2470" class="blob-code js-file-line">                        undbzValue = (<span class="pl-k">float</span>) <span class="pl-c1">exp</span>(<span class="pl-c1">0.1</span>*<span class="pl-c1">log</span>(<span class="pl-c1">10</span>)*dbzValue);</td>
      </tr>
      <tr>
        <td id="L2471" class="blob-num js-line-number" data-line-number="2471"></td>
        <td id="LC2471" class="blob-code js-file-line">                        <span class="pl-c">// sum the undbz in this layer</span></td>
      </tr>
      <tr>
        <td id="L2472" class="blob-num js-line-number" data-line-number="2472"></td>
        <td id="LC2472" class="blob-code js-file-line">                        undbzSum += undbzValue;</td>
      </tr>
      <tr>
        <td id="L2473" class="blob-num js-line-number" data-line-number="2473"></td>
        <td id="LC2473" class="blob-code js-file-line">                        <span class="pl-c">// copy the observed vrad value at this [azimuth, elevation] </span></td>
      </tr>
      <tr>
        <td id="L2474" class="blob-num js-line-number" data-line-number="2474"></td>
        <td id="LC2474" class="blob-code js-file-line">                        yObs[iPointIncluded] = points[iPointLayer * nColsPoints + vradValueCol];</td>
      </tr>
      <tr>
        <td id="L2475" class="blob-num js-line-number" data-line-number="2475"></td>
        <td id="LC2475" class="blob-code js-file-line">                        <span class="pl-c">// pre-allocate the fitted vrad value at this [azimuth,elevation]</span></td>
      </tr>
      <tr>
        <td id="L2476" class="blob-num js-line-number" data-line-number="2476"></td>
        <td id="LC2476" class="blob-code js-file-line">                        yFitted[iPointIncluded] = <span class="pl-c1">0</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L2477" class="blob-num js-line-number" data-line-number="2477"></td>
        <td id="LC2477" class="blob-code js-file-line">                        <span class="pl-c">// keep a record of which index was just included</span></td>
      </tr>
      <tr>
        <td id="L2478" class="blob-num js-line-number" data-line-number="2478"></td>
        <td id="LC2478" class="blob-code js-file-line">                        includedIndex[iPointIncluded] = iPointLayer;</td>
      </tr>
      <tr>
        <td id="L2479" class="blob-num js-line-number" data-line-number="2479"></td>
        <td id="LC2479" class="blob-code js-file-line">                        </td>
      </tr>
      <tr>
        <td id="L2480" class="blob-num js-line-number" data-line-number="2480"></td>
        <td id="LC2480" class="blob-code js-file-line">                        <span class="pl-c">// raise the counter</span></td>
      </tr>
      <tr>
        <td id="L2481" class="blob-num js-line-number" data-line-number="2481"></td>
        <td id="LC2481" class="blob-code js-file-line">                        iPointIncluded += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2482" class="blob-num js-line-number" data-line-number="2482"></td>
        <td id="LC2482" class="blob-code js-file-line">                        </td>
      </tr>
      <tr>
        <td id="L2483" class="blob-num js-line-number" data-line-number="2483"></td>
        <td id="LC2483" class="blob-code js-file-line">                    }</td>
      </tr>
      <tr>
        <td id="L2484" class="blob-num js-line-number" data-line-number="2484"></td>
        <td id="LC2484" class="blob-code js-file-line">                } <span class="pl-c">// endfor (iPointLayer = 0; iPointLayer &lt; nPointsLayer; iPointLayer++) {</span></td>
      </tr>
      <tr>
        <td id="L2485" class="blob-num js-line-number" data-line-number="2485"></td>
        <td id="LC2485" class="blob-code js-file-line">                nPointsIncluded = iPointIncluded;</td>
      </tr>
      <tr>
        <td id="L2486" class="blob-num js-line-number" data-line-number="2486"></td>
        <td id="LC2486" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2487" class="blob-num js-line-number" data-line-number="2487"></td>
        <td id="LC2487" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2488" class="blob-num js-line-number" data-line-number="2488"></td>
        <td id="LC2488" class="blob-code js-file-line">                <span class="pl-k">if</span> (nPointsIncluded &gt; nPointsIncludedMin) {   </td>
      </tr>
      <tr>
        <td id="L2489" class="blob-num js-line-number" data-line-number="2489"></td>
        <td id="LC2489" class="blob-code js-file-line">                    <span class="pl-c">// when there are enough valid points, convert undbzAvg back to dB-scale</span></td>
      </tr>
      <tr>
        <td id="L2490" class="blob-num js-line-number" data-line-number="2490"></td>
        <td id="LC2490" class="blob-code js-file-line">                    undbzAvg = (<span class="pl-k">float</span>) (undbzSum/nPointsIncluded);</td>
      </tr>
      <tr>
        <td id="L2491" class="blob-num js-line-number" data-line-number="2491"></td>
        <td id="LC2491" class="blob-code js-file-line">                    dbzAvg = (<span class="pl-c1">10</span>*<span class="pl-c1">log</span>(undbzAvg))/<span class="pl-c1">log</span>(<span class="pl-c1">10</span>);</td>
      </tr>
      <tr>
        <td id="L2492" class="blob-num js-line-number" data-line-number="2492"></td>
        <td id="LC2492" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L2493" class="blob-num js-line-number" data-line-number="2493"></td>
        <td id="LC2493" class="blob-code js-file-line">                <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2494" class="blob-num js-line-number" data-line-number="2494"></td>
        <td id="LC2494" class="blob-code js-file-line">                    undbzAvg = NAN;</td>
      </tr>
      <tr>
        <td id="L2495" class="blob-num js-line-number" data-line-number="2495"></td>
        <td id="LC2495" class="blob-code js-file-line">                    dbzAvg = NAN;</td>
      </tr>
      <tr>
        <td id="L2496" class="blob-num js-line-number" data-line-number="2496"></td>
        <td id="LC2496" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L2497" class="blob-num js-line-number" data-line-number="2497"></td>
        <td id="LC2497" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2498" class="blob-num js-line-number" data-line-number="2498"></td>
        <td id="LC2498" class="blob-code js-file-line">                <span class="pl-c">// convert from Z (not dBZ) in units of mm^6/m^3 to </span></td>
      </tr>
      <tr>
        <td id="L2499" class="blob-num js-line-number" data-line-number="2499"></td>
        <td id="LC2499" class="blob-code js-file-line">                <span class="pl-c">// reflectivity eta in units of cm^2/km^3</span></td>
      </tr>
      <tr>
        <td id="L2500" class="blob-num js-line-number" data-line-number="2500"></td>
        <td id="LC2500" class="blob-code js-file-line">                reflectivity = dbzFactor * undbzAvg;</td>
      </tr>
      <tr>
        <td id="L2501" class="blob-num js-line-number" data-line-number="2501"></td>
        <td id="LC2501" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2502" class="blob-num js-line-number" data-line-number="2502"></td>
        <td id="LC2502" class="blob-code js-file-line">                <span class="pl-k">if</span> (iProfileType == <span class="pl-c1">1</span>) {</td>
      </tr>
      <tr>
        <td id="L2503" class="blob-num js-line-number" data-line-number="2503"></td>
        <td id="LC2503" class="blob-code js-file-line">                    <span class="pl-c">// calculate bird density in number of birds/km^3 by</span></td>
      </tr>
      <tr>
        <td id="L2504" class="blob-num js-line-number" data-line-number="2504"></td>
        <td id="LC2504" class="blob-code js-file-line">                    <span class="pl-c">// dividing the reflectivity by the (assumed) cross section</span></td>
      </tr>
      <tr>
        <td id="L2505" class="blob-num js-line-number" data-line-number="2505"></td>
        <td id="LC2505" class="blob-code js-file-line">                    <span class="pl-c">// of one bird</span></td>
      </tr>
      <tr>
        <td id="L2506" class="blob-num js-line-number" data-line-number="2506"></td>
        <td id="LC2506" class="blob-code js-file-line">                    birdDensity = reflectivity / birdRadarCrossSection;</td>
      </tr>
      <tr>
        <td id="L2507" class="blob-num js-line-number" data-line-number="2507"></td>
        <td id="LC2507" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L2508" class="blob-num js-line-number" data-line-number="2508"></td>
        <td id="LC2508" class="blob-code js-file-line">                <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2509" class="blob-num js-line-number" data-line-number="2509"></td>
        <td id="LC2509" class="blob-code js-file-line">                    birdDensity = NAN;</td>
      </tr>
      <tr>
        <td id="L2510" class="blob-num js-line-number" data-line-number="2510"></td>
        <td id="LC2510" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L2511" class="blob-num js-line-number" data-line-number="2511"></td>
        <td id="LC2511" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2512" class="blob-num js-line-number" data-line-number="2512"></td>
        <td id="LC2512" class="blob-code js-file-line">                <span class="pl-c">// check if there are directions that have almost no observations</span></td>
      </tr>
      <tr>
        <td id="L2513" class="blob-num js-line-number" data-line-number="2513"></td>
        <td id="LC2513" class="blob-code js-file-line">                <span class="pl-c">// (as this makes the svdfit result really uncertain)  </span></td>
      </tr>
      <tr>
        <td id="L2514" class="blob-num js-line-number" data-line-number="2514"></td>
        <td id="LC2514" class="blob-code js-file-line">                hasGap = <span class="pl-c1">hasAzimuthGap</span>(&amp;pointsSelection[<span class="pl-c1">0</span>], nPointsIncluded);</td>
      </tr>
      <tr>
        <td id="L2515" class="blob-num js-line-number" data-line-number="2515"></td>
        <td id="LC2515" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2516" class="blob-num js-line-number" data-line-number="2516"></td>
        <td id="LC2516" class="blob-code js-file-line">                <span class="pl-k">if</span> (fitVrad == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L2517" class="blob-num js-line-number" data-line-number="2517"></td>
        <td id="LC2517" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2518" class="blob-num js-line-number" data-line-number="2518"></td>
        <td id="LC2518" class="blob-code js-file-line">                    <span class="pl-k">if</span> (hasGap==<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L2519" class="blob-num js-line-number" data-line-number="2519"></td>
        <td id="LC2519" class="blob-code js-file-line">                        </td>
      </tr>
      <tr>
        <td id="L2520" class="blob-num js-line-number" data-line-number="2520"></td>
        <td id="LC2520" class="blob-code js-file-line">                        <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2521" class="blob-num js-line-number" data-line-number="2521"></td>
        <td id="LC2521" class="blob-code js-file-line">                        <span class="pl-c">//                       do the svdfit                           //</span></td>
      </tr>
      <tr>
        <td id="L2522" class="blob-num js-line-number" data-line-number="2522"></td>
        <td id="LC2522" class="blob-code js-file-line">                        <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2523" class="blob-num js-line-number" data-line-number="2523"></td>
        <td id="LC2523" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2524" class="blob-num js-line-number" data-line-number="2524"></td>
        <td id="LC2524" class="blob-code js-file-line">                        chisq = <span class="pl-c1">svdfit</span>(&amp;pointsSelection[<span class="pl-c1">0</span>], nDims, &amp;yObs[<span class="pl-c1">0</span>], &amp;yFitted[<span class="pl-c1">0</span>], </td>
      </tr>
      <tr>
        <td id="L2525" class="blob-num js-line-number" data-line-number="2525"></td>
        <td id="LC2525" class="blob-code js-file-line">                                nPointsIncluded, &amp;parameterVector[<span class="pl-c1">0</span>], &amp;avar[<span class="pl-c1">0</span>], nParsFitted);</td>
      </tr>
      <tr>
        <td id="L2526" class="blob-num js-line-number" data-line-number="2526"></td>
        <td id="LC2526" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2527" class="blob-num js-line-number" data-line-number="2527"></td>
        <td id="LC2527" class="blob-code js-file-line">                        <span class="pl-k">if</span> (chisq &lt; chisqMin) {</td>
      </tr>
      <tr>
        <td id="L2528" class="blob-num js-line-number" data-line-number="2528"></td>
        <td id="LC2528" class="blob-code js-file-line">                            <span class="pl-c">// the fit was not good enough, reset parameter vector array </span></td>
      </tr>
      <tr>
        <td id="L2529" class="blob-num js-line-number" data-line-number="2529"></td>
        <td id="LC2529" class="blob-code js-file-line">                            <span class="pl-c">// elements to NAN and continue with the next layer</span></td>
      </tr>
      <tr>
        <td id="L2530" class="blob-num js-line-number" data-line-number="2530"></td>
        <td id="LC2530" class="blob-code js-file-line">                            parameterVector[<span class="pl-c1">0</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2531" class="blob-num js-line-number" data-line-number="2531"></td>
        <td id="LC2531" class="blob-code js-file-line">                            parameterVector[<span class="pl-c1">1</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2532" class="blob-num js-line-number" data-line-number="2532"></td>
        <td id="LC2532" class="blob-code js-file-line">                            parameterVector[<span class="pl-c1">2</span>] = NAN;</td>
      </tr>
      <tr>
        <td id="L2533" class="blob-num js-line-number" data-line-number="2533"></td>
        <td id="LC2533" class="blob-code js-file-line">                            <span class="pl-k">continue</span>;</td>
      </tr>
      <tr>
        <td id="L2534" class="blob-num js-line-number" data-line-number="2534"></td>
        <td id="LC2534" class="blob-code js-file-line">                        } </td>
      </tr>
      <tr>
        <td id="L2535" class="blob-num js-line-number" data-line-number="2535"></td>
        <td id="LC2535" class="blob-code js-file-line">                        <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2536" class="blob-num js-line-number" data-line-number="2536"></td>
        <td id="LC2536" class="blob-code js-file-line">                            </td>
      </tr>
      <tr>
        <td id="L2537" class="blob-num js-line-number" data-line-number="2537"></td>
        <td id="LC2537" class="blob-code js-file-line">                            chi = <span class="pl-c1">sqrt</span>(chisq);</td>
      </tr>
      <tr>
        <td id="L2538" class="blob-num js-line-number" data-line-number="2538"></td>
        <td id="LC2538" class="blob-code js-file-line">                            hSpeed = <span class="pl-c1">sqrt</span>(<span class="pl-c1">pow</span>(parameterVector[<span class="pl-c1">0</span>],<span class="pl-c1">2</span>) + <span class="pl-c1">pow</span>(parameterVector[<span class="pl-c1">1</span>],<span class="pl-c1">2</span>));</td>
      </tr>
      <tr>
        <td id="L2539" class="blob-num js-line-number" data-line-number="2539"></td>
        <td id="LC2539" class="blob-code js-file-line">                            hDir = (<span class="pl-c1">atan2</span>(parameterVector[<span class="pl-c1">0</span>],parameterVector[<span class="pl-c1">1</span>])*RAD2DEG);</td>
      </tr>
      <tr>
        <td id="L2540" class="blob-num js-line-number" data-line-number="2540"></td>
        <td id="LC2540" class="blob-code js-file-line">                            </td>
      </tr>
      <tr>
        <td id="L2541" class="blob-num js-line-number" data-line-number="2541"></td>
        <td id="LC2541" class="blob-code js-file-line">                            <span class="pl-k">if</span> (hDir &lt; <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L2542" class="blob-num js-line-number" data-line-number="2542"></td>
        <td id="LC2542" class="blob-code js-file-line">                                hDir += <span class="pl-c1">360</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L2543" class="blob-num js-line-number" data-line-number="2543"></td>
        <td id="LC2543" class="blob-code js-file-line">                            }</td>
      </tr>
      <tr>
        <td id="L2544" class="blob-num js-line-number" data-line-number="2544"></td>
        <td id="LC2544" class="blob-code js-file-line">                            </td>
      </tr>
      <tr>
        <td id="L2545" class="blob-num js-line-number" data-line-number="2545"></td>
        <td id="LC2545" class="blob-code js-file-line">                        }</td>
      </tr>
      <tr>
        <td id="L2546" class="blob-num js-line-number" data-line-number="2546"></td>
        <td id="LC2546" class="blob-code js-file-line">                    } <span class="pl-c">// endif (hasGap == FALSE)</span></td>
      </tr>
      <tr>
        <td id="L2547" class="blob-num js-line-number" data-line-number="2547"></td>
        <td id="LC2547" class="blob-code js-file-line">                    </td>
      </tr>
      <tr>
        <td id="L2548" class="blob-num js-line-number" data-line-number="2548"></td>
        <td id="LC2548" class="blob-code js-file-line">                    <span class="pl-c">// if the fitted vrad value is more than &#39;absVDifMax&#39; away from the corresponding</span></td>
      </tr>
      <tr>
        <td id="L2549" class="blob-num js-line-number" data-line-number="2549"></td>
        <td id="LC2549" class="blob-code js-file-line">                    <span class="pl-c">// observed vrad value, set the gate&#39;s flagPositionVDifMax bit flag to 1, excluding the </span></td>
      </tr>
      <tr>
        <td id="L2550" class="blob-num js-line-number" data-line-number="2550"></td>
        <td id="LC2550" class="blob-code js-file-line">                    <span class="pl-c">// gate in the second svdfit iteration</span></td>
      </tr>
      <tr>
        <td id="L2551" class="blob-num js-line-number" data-line-number="2551"></td>
        <td id="LC2551" class="blob-code js-file-line">                    <span class="pl-c1">updateFlagFieldsInPointsArray</span>(&amp;yObs[<span class="pl-c1">0</span>], &amp;yFitted[<span class="pl-c1">0</span>], &amp;includedIndex[<span class="pl-c1">0</span>], nPointsIncluded,</td>
      </tr>
      <tr>
        <td id="L2552" class="blob-num js-line-number" data-line-number="2552"></td>
        <td id="LC2552" class="blob-code js-file-line">                                                  &amp;points[<span class="pl-c1">0</span>]);</td>
      </tr>
      <tr>
        <td id="L2553" class="blob-num js-line-number" data-line-number="2553"></td>
        <td id="LC2553" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2554" class="blob-num js-line-number" data-line-number="2554"></td>
        <td id="LC2554" class="blob-code js-file-line">                }; <span class="pl-c">// endif (fitVrad == TRUE)</span></td>
      </tr>
      <tr>
        <td id="L2555" class="blob-num js-line-number" data-line-number="2555"></td>
        <td id="LC2555" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2556" class="blob-num js-line-number" data-line-number="2556"></td>
        <td id="LC2556" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">0</span>] = iLayer * layerThickness;</td>
      </tr>
      <tr>
        <td id="L2557" class="blob-num js-line-number" data-line-number="2557"></td>
        <td id="LC2557" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">1</span>] = (iLayer + <span class="pl-c1">1</span>) * layerThickness;</td>
      </tr>
      <tr>
        <td id="L2558" class="blob-num js-line-number" data-line-number="2558"></td>
        <td id="LC2558" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">2</span>] = parameterVector[<span class="pl-c1">0</span>];</td>
      </tr>
      <tr>
        <td id="L2559" class="blob-num js-line-number" data-line-number="2559"></td>
        <td id="LC2559" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">3</span>] = parameterVector[<span class="pl-c1">1</span>];</td>
      </tr>
      <tr>
        <td id="L2560" class="blob-num js-line-number" data-line-number="2560"></td>
        <td id="LC2560" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">4</span>] = parameterVector[<span class="pl-c1">2</span>];</td>
      </tr>
      <tr>
        <td id="L2561" class="blob-num js-line-number" data-line-number="2561"></td>
        <td id="LC2561" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">5</span>] = hSpeed;</td>
      </tr>
      <tr>
        <td id="L2562" class="blob-num js-line-number" data-line-number="2562"></td>
        <td id="LC2562" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">6</span>] = hDir;</td>
      </tr>
      <tr>
        <td id="L2563" class="blob-num js-line-number" data-line-number="2563"></td>
        <td id="LC2563" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">7</span>] = chi;</td>
      </tr>
      <tr>
        <td id="L2564" class="blob-num js-line-number" data-line-number="2564"></td>
        <td id="LC2564" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">8</span>] = (<span class="pl-k">float</span>) hasGap;</td>
      </tr>
      <tr>
        <td id="L2565" class="blob-num js-line-number" data-line-number="2565"></td>
        <td id="LC2565" class="blob-code js-file-line">                profile[iLayer*nColsProfile +  <span class="pl-c1">9</span>] = dbzAvg;</td>
      </tr>
      <tr>
        <td id="L2566" class="blob-num js-line-number" data-line-number="2566"></td>
        <td id="LC2566" class="blob-code js-file-line">                profile[iLayer*nColsProfile + <span class="pl-c1">10</span>] = (<span class="pl-k">float</span>) nPointsIncluded;</td>
      </tr>
      <tr>
        <td id="L2567" class="blob-num js-line-number" data-line-number="2567"></td>
        <td id="LC2567" class="blob-code js-file-line">                profile[iLayer*nColsProfile + <span class="pl-c1">11</span>] = reflectivity;</td>
      </tr>
      <tr>
        <td id="L2568" class="blob-num js-line-number" data-line-number="2568"></td>
        <td id="LC2568" class="blob-code js-file-line">                profile[iLayer*nColsProfile + <span class="pl-c1">12</span>] = birdDensity;</td>
      </tr>
      <tr>
        <td id="L2569" class="blob-num js-line-number" data-line-number="2569"></td>
        <td id="LC2569" class="blob-code js-file-line">                </td>
      </tr>
      <tr>
        <td id="L2570" class="blob-num js-line-number" data-line-number="2570"></td>
        <td id="LC2570" class="blob-code js-file-line">                <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) yObs);</td>
      </tr>
      <tr>
        <td id="L2571" class="blob-num js-line-number" data-line-number="2571"></td>
        <td id="LC2571" class="blob-code js-file-line">                <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) yFitted);</td>
      </tr>
      <tr>
        <td id="L2572" class="blob-num js-line-number" data-line-number="2572"></td>
        <td id="LC2572" class="blob-code js-file-line">                <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) pointsSelection);</td>
      </tr>
      <tr>
        <td id="L2573" class="blob-num js-line-number" data-line-number="2573"></td>
        <td id="LC2573" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2574" class="blob-num js-line-number" data-line-number="2574"></td>
        <td id="LC2574" class="blob-code js-file-line">            } <span class="pl-c">// endfor (iPass = 0; iPass &lt; nPasses; iPass++)</span></td>
      </tr>
      <tr>
        <td id="L2575" class="blob-num js-line-number" data-line-number="2575"></td>
        <td id="LC2575" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2576" class="blob-num js-line-number" data-line-number="2576"></td>
        <td id="LC2576" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2577" class="blob-num js-line-number" data-line-number="2577"></td>
        <td id="LC2577" class="blob-code js-file-line">            <span class="pl-c">// You need some of the results of iProfileType == 3 in order </span></td>
      </tr>
      <tr>
        <td id="L2578" class="blob-num js-line-number" data-line-number="2578"></td>
        <td id="LC2578" class="blob-code js-file-line">            <span class="pl-c">// to calculate iProfileType == 1</span></td>
      </tr>
      <tr>
        <td id="L2579" class="blob-num js-line-number" data-line-number="2579"></td>
        <td id="LC2579" class="blob-code js-file-line">            <span class="pl-k">if</span> (iProfileType == <span class="pl-c1">3</span>) {</td>
      </tr>
      <tr>
        <td id="L2580" class="blob-num js-line-number" data-line-number="2580"></td>
        <td id="LC2580" class="blob-code js-file-line">                <span class="pl-k">if</span> (chi &lt; stdDevMinBird) {</td>
      </tr>
      <tr>
        <td id="L2581" class="blob-num js-line-number" data-line-number="2581"></td>
        <td id="LC2581" class="blob-code js-file-line">                    scatterersAreNotBirds[iLayer] = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L2582" class="blob-num js-line-number" data-line-number="2582"></td>
        <td id="LC2582" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L2583" class="blob-num js-line-number" data-line-number="2583"></td>
        <td id="LC2583" class="blob-code js-file-line">                <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L2584" class="blob-num js-line-number" data-line-number="2584"></td>
        <td id="LC2584" class="blob-code js-file-line">                    scatterersAreNotBirds[iLayer] = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L2585" class="blob-num js-line-number" data-line-number="2585"></td>
        <td id="LC2585" class="blob-code js-file-line">                }</td>
      </tr>
      <tr>
        <td id="L2586" class="blob-num js-line-number" data-line-number="2586"></td>
        <td id="LC2586" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L2587" class="blob-num js-line-number" data-line-number="2587"></td>
        <td id="LC2587" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2588" class="blob-num js-line-number" data-line-number="2588"></td>
        <td id="LC2588" class="blob-code js-file-line">        } <span class="pl-c">// endfor (iLayer = 0; iLayer &lt; nLayers; iLayer++)</span></td>
      </tr>
      <tr>
        <td id="L2589" class="blob-num js-line-number" data-line-number="2589"></td>
        <td id="LC2589" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2590" class="blob-num js-line-number" data-line-number="2590"></td>
        <td id="LC2590" class="blob-code js-file-line">        <span class="pl-k">if</span> (printProfileVar ==  <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L2591" class="blob-num js-line-number" data-line-number="2591"></td>
        <td id="LC2591" class="blob-code js-file-line">            <span class="pl-c1">printProfile</span>();</td>
      </tr>
      <tr>
        <td id="L2592" class="blob-num js-line-number" data-line-number="2592"></td>
        <td id="LC2592" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2593" class="blob-num js-line-number" data-line-number="2593"></td>
        <td id="LC2593" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2594" class="blob-num js-line-number" data-line-number="2594"></td>
        <td id="LC2594" class="blob-code js-file-line">    } <span class="pl-c">// endfor (iProfileType = nProfileTypes; iProfileType &gt; 0; iProfileType--)</span></td>
      </tr>
      <tr>
        <td id="L2595" class="blob-num js-line-number" data-line-number="2595"></td>
        <td id="LC2595" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2596" class="blob-num js-line-number" data-line-number="2596"></td>
        <td id="LC2596" class="blob-code js-file-line">} <span class="pl-c">// vol2birdCalcProfiles</span></td>
      </tr>
      <tr>
        <td id="L2597" class="blob-num js-line-number" data-line-number="2597"></td>
        <td id="LC2597" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2598" class="blob-num js-line-number" data-line-number="2598"></td>
        <td id="LC2598" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2599" class="blob-num js-line-number" data-line-number="2599"></td>
        <td id="LC2599" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2600" class="blob-num js-line-number" data-line-number="2600"></td>
        <td id="LC2600" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2601" class="blob-num js-line-number" data-line-number="2601"></td>
        <td id="LC2601" class="blob-code js-file-line"><span class="pl-k">void</span> <span class="pl-en">printImageInt</span>(<span class="pl-k">const</span> SCANMETA* meta, <span class="pl-k">const</span> <span class="pl-k">int</span>* imageInt) {</td>
      </tr>
      <tr>
        <td id="L2602" class="blob-num js-line-number" data-line-number="2602"></td>
        <td id="LC2602" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2603" class="blob-num js-line-number" data-line-number="2603"></td>
        <td id="LC2603" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2604" class="blob-num js-line-number" data-line-number="2604"></td>
        <td id="LC2604" class="blob-code js-file-line">    <span class="pl-k">int</span> nRang = meta-&gt;nRang;</td>
      </tr>
      <tr>
        <td id="L2605" class="blob-num js-line-number" data-line-number="2605"></td>
        <td id="LC2605" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzim = meta-&gt;nAzim;</td>
      </tr>
      <tr>
        <td id="L2606" class="blob-num js-line-number" data-line-number="2606"></td>
        <td id="LC2606" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L2607" class="blob-num js-line-number" data-line-number="2607"></td>
        <td id="LC2607" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzim;</td>
      </tr>
      <tr>
        <td id="L2608" class="blob-num js-line-number" data-line-number="2608"></td>
        <td id="LC2608" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L2609" class="blob-num js-line-number" data-line-number="2609"></td>
        <td id="LC2609" class="blob-code js-file-line">    <span class="pl-k">int</span> needsSignChar;</td>
      </tr>
      <tr>
        <td id="L2610" class="blob-num js-line-number" data-line-number="2610"></td>
        <td id="LC2610" class="blob-code js-file-line">    <span class="pl-k">int</span> maxValue;</td>
      </tr>
      <tr>
        <td id="L2611" class="blob-num js-line-number" data-line-number="2611"></td>
        <td id="LC2611" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2612" class="blob-num js-line-number" data-line-number="2612"></td>
        <td id="LC2612" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2613" class="blob-num js-line-number" data-line-number="2613"></td>
        <td id="LC2613" class="blob-code js-file-line">    <span class="pl-k">int</span> thisValue;</td>
      </tr>
      <tr>
        <td id="L2614" class="blob-num js-line-number" data-line-number="2614"></td>
        <td id="LC2614" class="blob-code js-file-line">    <span class="pl-k">int</span> nChars;</td>
      </tr>
      <tr>
        <td id="L2615" class="blob-num js-line-number" data-line-number="2615"></td>
        <td id="LC2615" class="blob-code js-file-line">    <span class="pl-k">char</span>* formatStr;</td>
      </tr>
      <tr>
        <td id="L2616" class="blob-num js-line-number" data-line-number="2616"></td>
        <td id="LC2616" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2617" class="blob-num js-line-number" data-line-number="2617"></td>
        <td id="LC2617" class="blob-code js-file-line">    iGlobal = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2618" class="blob-num js-line-number" data-line-number="2618"></td>
        <td id="LC2618" class="blob-code js-file-line">    maxValue = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2619" class="blob-num js-line-number" data-line-number="2619"></td>
        <td id="LC2619" class="blob-code js-file-line">    needsSignChar = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L2620" class="blob-num js-line-number" data-line-number="2620"></td>
        <td id="LC2620" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2621" class="blob-num js-line-number" data-line-number="2621"></td>
        <td id="LC2621" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>elevAngle = <span class="pl-c1">%f</span> degrees<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,meta-&gt;elev);</td>
      </tr>
      <tr>
        <td id="L2622" class="blob-num js-line-number" data-line-number="2622"></td>
        <td id="LC2622" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2623" class="blob-num js-line-number" data-line-number="2623"></td>
        <td id="LC2623" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2624" class="blob-num js-line-number" data-line-number="2624"></td>
        <td id="LC2624" class="blob-code js-file-line">    <span class="pl-c">// first, determine how many characters are needed to print array &#39;imageInt&#39;</span></td>
      </tr>
      <tr>
        <td id="L2625" class="blob-num js-line-number" data-line-number="2625"></td>
        <td id="LC2625" class="blob-code js-file-line">    <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) { </td>
      </tr>
      <tr>
        <td id="L2626" class="blob-num js-line-number" data-line-number="2626"></td>
        <td id="LC2626" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2627" class="blob-num js-line-number" data-line-number="2627"></td>
        <td id="LC2627" class="blob-code js-file-line">        <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L2628" class="blob-num js-line-number" data-line-number="2628"></td>
        <td id="LC2628" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2629" class="blob-num js-line-number" data-line-number="2629"></td>
        <td id="LC2629" class="blob-code js-file-line">            thisValue = (<span class="pl-k">int</span>) <span class="pl-c1">XABS</span>(imageInt[iGlobal]);</td>
      </tr>
      <tr>
        <td id="L2630" class="blob-num js-line-number" data-line-number="2630"></td>
        <td id="LC2630" class="blob-code js-file-line">            <span class="pl-k">if</span> (imageInt[iGlobal] &lt; <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L2631" class="blob-num js-line-number" data-line-number="2631"></td>
        <td id="LC2631" class="blob-code js-file-line">                needsSignChar = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L2632" class="blob-num js-line-number" data-line-number="2632"></td>
        <td id="LC2632" class="blob-code js-file-line">            }</td>
      </tr>
      <tr>
        <td id="L2633" class="blob-num js-line-number" data-line-number="2633"></td>
        <td id="LC2633" class="blob-code js-file-line">            <span class="pl-k">if</span> (thisValue &gt; maxValue) {</td>
      </tr>
      <tr>
        <td id="L2634" class="blob-num js-line-number" data-line-number="2634"></td>
        <td id="LC2634" class="blob-code js-file-line">                maxValue = thisValue;</td>
      </tr>
      <tr>
        <td id="L2635" class="blob-num js-line-number" data-line-number="2635"></td>
        <td id="LC2635" class="blob-code js-file-line">            } ;</td>
      </tr>
      <tr>
        <td id="L2636" class="blob-num js-line-number" data-line-number="2636"></td>
        <td id="LC2636" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2637" class="blob-num js-line-number" data-line-number="2637"></td>
        <td id="LC2637" class="blob-code js-file-line">            iGlobal += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2638" class="blob-num js-line-number" data-line-number="2638"></td>
        <td id="LC2638" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2639" class="blob-num js-line-number" data-line-number="2639"></td>
        <td id="LC2639" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2640" class="blob-num js-line-number" data-line-number="2640"></td>
        <td id="LC2640" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2641" class="blob-num js-line-number" data-line-number="2641"></td>
        <td id="LC2641" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2642" class="blob-num js-line-number" data-line-number="2642"></td>
        <td id="LC2642" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2643" class="blob-num js-line-number" data-line-number="2643"></td>
        <td id="LC2643" class="blob-code js-file-line">    nChars = (<span class="pl-k">int</span>) <span class="pl-c1">ceil</span>(<span class="pl-c1">log</span>(maxValue + <span class="pl-c1">1</span>)/<span class="pl-c1">log</span>(<span class="pl-c1">10</span>));</td>
      </tr>
      <tr>
        <td id="L2644" class="blob-num js-line-number" data-line-number="2644"></td>
        <td id="LC2644" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2645" class="blob-num js-line-number" data-line-number="2645"></td>
        <td id="LC2645" class="blob-code js-file-line">    <span class="pl-k">if</span> (needsSignChar) {</td>
      </tr>
      <tr>
        <td id="L2646" class="blob-num js-line-number" data-line-number="2646"></td>
        <td id="LC2646" class="blob-code js-file-line">        nChars += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2647" class="blob-num js-line-number" data-line-number="2647"></td>
        <td id="LC2647" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2648" class="blob-num js-line-number" data-line-number="2648"></td>
        <td id="LC2648" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2649" class="blob-num js-line-number" data-line-number="2649"></td>
        <td id="LC2649" class="blob-code js-file-line">    <span class="pl-k">switch</span> (nChars) {</td>
      </tr>
      <tr>
        <td id="L2650" class="blob-num js-line-number" data-line-number="2650"></td>
        <td id="LC2650" class="blob-code js-file-line">        <span class="pl-k">case</span> <span class="pl-c1">0</span> :</td>
      </tr>
      <tr>
        <td id="L2651" class="blob-num js-line-number" data-line-number="2651"></td>
        <td id="LC2651" class="blob-code js-file-line">            formatStr = <span class="pl-s"><span class="pl-pds">&quot;</span> <span class="pl-c1">%1d</span><span class="pl-pds">&quot;</span></span>;</td>
      </tr>
      <tr>
        <td id="L2652" class="blob-num js-line-number" data-line-number="2652"></td>
        <td id="LC2652" class="blob-code js-file-line">            <span class="pl-k">break</span>; </td>
      </tr>
      <tr>
        <td id="L2653" class="blob-num js-line-number" data-line-number="2653"></td>
        <td id="LC2653" class="blob-code js-file-line">        <span class="pl-k">case</span> <span class="pl-c1">1</span> :</td>
      </tr>
      <tr>
        <td id="L2654" class="blob-num js-line-number" data-line-number="2654"></td>
        <td id="LC2654" class="blob-code js-file-line">            formatStr = <span class="pl-s"><span class="pl-pds">&quot;</span> <span class="pl-c1">%1d</span><span class="pl-pds">&quot;</span></span>;</td>
      </tr>
      <tr>
        <td id="L2655" class="blob-num js-line-number" data-line-number="2655"></td>
        <td id="LC2655" class="blob-code js-file-line">            <span class="pl-k">break</span>; </td>
      </tr>
      <tr>
        <td id="L2656" class="blob-num js-line-number" data-line-number="2656"></td>
        <td id="LC2656" class="blob-code js-file-line">        <span class="pl-k">case</span> <span class="pl-c1">2</span> :</td>
      </tr>
      <tr>
        <td id="L2657" class="blob-num js-line-number" data-line-number="2657"></td>
        <td id="LC2657" class="blob-code js-file-line">            formatStr = <span class="pl-s"><span class="pl-pds">&quot;</span> <span class="pl-c1">%2d</span><span class="pl-pds">&quot;</span></span>; </td>
      </tr>
      <tr>
        <td id="L2658" class="blob-num js-line-number" data-line-number="2658"></td>
        <td id="LC2658" class="blob-code js-file-line">            <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L2659" class="blob-num js-line-number" data-line-number="2659"></td>
        <td id="LC2659" class="blob-code js-file-line">        <span class="pl-k">case</span> <span class="pl-c1">3</span> :</td>
      </tr>
      <tr>
        <td id="L2660" class="blob-num js-line-number" data-line-number="2660"></td>
        <td id="LC2660" class="blob-code js-file-line">            formatStr = <span class="pl-s"><span class="pl-pds">&quot;</span> <span class="pl-c1">%3d</span><span class="pl-pds">&quot;</span></span>; </td>
      </tr>
      <tr>
        <td id="L2661" class="blob-num js-line-number" data-line-number="2661"></td>
        <td id="LC2661" class="blob-code js-file-line">            <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L2662" class="blob-num js-line-number" data-line-number="2662"></td>
        <td id="LC2662" class="blob-code js-file-line">        <span class="pl-k">case</span> <span class="pl-c1">4</span> :</td>
      </tr>
      <tr>
        <td id="L2663" class="blob-num js-line-number" data-line-number="2663"></td>
        <td id="LC2663" class="blob-code js-file-line">            formatStr = <span class="pl-s"><span class="pl-pds">&quot;</span> <span class="pl-c1">%4d</span><span class="pl-pds">&quot;</span></span>; </td>
      </tr>
      <tr>
        <td id="L2664" class="blob-num js-line-number" data-line-number="2664"></td>
        <td id="LC2664" class="blob-code js-file-line">            <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L2665" class="blob-num js-line-number" data-line-number="2665"></td>
        <td id="LC2665" class="blob-code js-file-line">        <span class="pl-k">default</span> :</td>
      </tr>
      <tr>
        <td id="L2666" class="blob-num js-line-number" data-line-number="2666"></td>
        <td id="LC2666" class="blob-code js-file-line">            formatStr = <span class="pl-s"><span class="pl-pds">&quot;</span> <span class="pl-c1">%8d</span><span class="pl-pds">&quot;</span></span>; </td>
      </tr>
      <tr>
        <td id="L2667" class="blob-num js-line-number" data-line-number="2667"></td>
        <td id="LC2667" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2668" class="blob-num js-line-number" data-line-number="2668"></td>
        <td id="LC2668" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2669" class="blob-num js-line-number" data-line-number="2669"></td>
        <td id="LC2669" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2670" class="blob-num js-line-number" data-line-number="2670"></td>
        <td id="LC2670" class="blob-code js-file-line">    iGlobal = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2671" class="blob-num js-line-number" data-line-number="2671"></td>
        <td id="LC2671" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2672" class="blob-num js-line-number" data-line-number="2672"></td>
        <td id="LC2672" class="blob-code js-file-line">    <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) { </td>
      </tr>
      <tr>
        <td id="L2673" class="blob-num js-line-number" data-line-number="2673"></td>
        <td id="LC2673" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2674" class="blob-num js-line-number" data-line-number="2674"></td>
        <td id="LC2674" class="blob-code js-file-line">        <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L2675" class="blob-num js-line-number" data-line-number="2675"></td>
        <td id="LC2675" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2676" class="blob-num js-line-number" data-line-number="2676"></td>
        <td id="LC2676" class="blob-code js-file-line">            iGlobal = iRang + iAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L2677" class="blob-num js-line-number" data-line-number="2677"></td>
        <td id="LC2677" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2678" class="blob-num js-line-number" data-line-number="2678"></td>
        <td id="LC2678" class="blob-code js-file-line">            thisValue = (<span class="pl-k">int</span>) imageInt[iGlobal];</td>
      </tr>
      <tr>
        <td id="L2679" class="blob-num js-line-number" data-line-number="2679"></td>
        <td id="LC2679" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2680" class="blob-num js-line-number" data-line-number="2680"></td>
        <td id="LC2680" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr,formatStr,thisValue);</td>
      </tr>
      <tr>
        <td id="L2681" class="blob-num js-line-number" data-line-number="2681"></td>
        <td id="LC2681" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2682" class="blob-num js-line-number" data-line-number="2682"></td>
        <td id="LC2682" class="blob-code js-file-line">            iGlobal += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2683" class="blob-num js-line-number" data-line-number="2683"></td>
        <td id="LC2683" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2684" class="blob-num js-line-number" data-line-number="2684"></td>
        <td id="LC2684" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2685" class="blob-num js-line-number" data-line-number="2685"></td>
        <td id="LC2685" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2686" class="blob-num js-line-number" data-line-number="2686"></td>
        <td id="LC2686" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2687" class="blob-num js-line-number" data-line-number="2687"></td>
        <td id="LC2687" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2688" class="blob-num js-line-number" data-line-number="2688"></td>
        <td id="LC2688" class="blob-code js-file-line">} <span class="pl-c">// printImageInt</span></td>
      </tr>
      <tr>
        <td id="L2689" class="blob-num js-line-number" data-line-number="2689"></td>
        <td id="LC2689" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2690" class="blob-num js-line-number" data-line-number="2690"></td>
        <td id="LC2690" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2691" class="blob-num js-line-number" data-line-number="2691"></td>
        <td id="LC2691" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2692" class="blob-num js-line-number" data-line-number="2692"></td>
        <td id="LC2692" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2693" class="blob-num js-line-number" data-line-number="2693"></td>
        <td id="LC2693" class="blob-code js-file-line"><span class="pl-k">void</span> <span class="pl-en">printImageUChar</span>(<span class="pl-k">const</span> SCANMETA* meta, <span class="pl-k">const</span> <span class="pl-k">unsigned</span> <span class="pl-k">char</span>* imageUChar) {</td>
      </tr>
      <tr>
        <td id="L2694" class="blob-num js-line-number" data-line-number="2694"></td>
        <td id="LC2694" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2695" class="blob-num js-line-number" data-line-number="2695"></td>
        <td id="LC2695" class="blob-code js-file-line">    <span class="pl-k">int</span> nAzim;</td>
      </tr>
      <tr>
        <td id="L2696" class="blob-num js-line-number" data-line-number="2696"></td>
        <td id="LC2696" class="blob-code js-file-line">    <span class="pl-k">int</span> iAzim;</td>
      </tr>
      <tr>
        <td id="L2697" class="blob-num js-line-number" data-line-number="2697"></td>
        <td id="LC2697" class="blob-code js-file-line">    <span class="pl-k">int</span> nRang;</td>
      </tr>
      <tr>
        <td id="L2698" class="blob-num js-line-number" data-line-number="2698"></td>
        <td id="LC2698" class="blob-code js-file-line">    <span class="pl-k">int</span> iRang;</td>
      </tr>
      <tr>
        <td id="L2699" class="blob-num js-line-number" data-line-number="2699"></td>
        <td id="LC2699" class="blob-code js-file-line">    <span class="pl-k">int</span> iGlobal;</td>
      </tr>
      <tr>
        <td id="L2700" class="blob-num js-line-number" data-line-number="2700"></td>
        <td id="LC2700" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2701" class="blob-num js-line-number" data-line-number="2701"></td>
        <td id="LC2701" class="blob-code js-file-line">    nAzim = meta-&gt;nAzim;</td>
      </tr>
      <tr>
        <td id="L2702" class="blob-num js-line-number" data-line-number="2702"></td>
        <td id="LC2702" class="blob-code js-file-line">    nRang = meta-&gt;nRang;</td>
      </tr>
      <tr>
        <td id="L2703" class="blob-num js-line-number" data-line-number="2703"></td>
        <td id="LC2703" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2704" class="blob-num js-line-number" data-line-number="2704"></td>
        <td id="LC2704" class="blob-code js-file-line">    <span class="pl-k">int</span>* imageInt = <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nRang * nAzim);</td>
      </tr>
      <tr>
        <td id="L2705" class="blob-num js-line-number" data-line-number="2705"></td>
        <td id="LC2705" class="blob-code js-file-line">    iGlobal = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2706" class="blob-num js-line-number" data-line-number="2706"></td>
        <td id="LC2706" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2707" class="blob-num js-line-number" data-line-number="2707"></td>
        <td id="LC2707" class="blob-code js-file-line">    <span class="pl-k">for</span> (iAzim = <span class="pl-c1">0</span>; iAzim &lt; nAzim; iAzim++) {</td>
      </tr>
      <tr>
        <td id="L2708" class="blob-num js-line-number" data-line-number="2708"></td>
        <td id="LC2708" class="blob-code js-file-line">        <span class="pl-k">for</span> (iRang = <span class="pl-c1">0</span>; iRang &lt; nRang; iRang++) {</td>
      </tr>
      <tr>
        <td id="L2709" class="blob-num js-line-number" data-line-number="2709"></td>
        <td id="LC2709" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2710" class="blob-num js-line-number" data-line-number="2710"></td>
        <td id="LC2710" class="blob-code js-file-line">            iGlobal = iRang + iAzim * nRang;</td>
      </tr>
      <tr>
        <td id="L2711" class="blob-num js-line-number" data-line-number="2711"></td>
        <td id="LC2711" class="blob-code js-file-line">            imageInt[iGlobal] = (<span class="pl-k">int</span>) imageUChar[iGlobal];</td>
      </tr>
      <tr>
        <td id="L2712" class="blob-num js-line-number" data-line-number="2712"></td>
        <td id="LC2712" class="blob-code js-file-line">            iGlobal += <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2713" class="blob-num js-line-number" data-line-number="2713"></td>
        <td id="LC2713" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2714" class="blob-num js-line-number" data-line-number="2714"></td>
        <td id="LC2714" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L2715" class="blob-num js-line-number" data-line-number="2715"></td>
        <td id="LC2715" class="blob-code js-file-line">    }     </td>
      </tr>
      <tr>
        <td id="L2716" class="blob-num js-line-number" data-line-number="2716"></td>
        <td id="LC2716" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2717" class="blob-num js-line-number" data-line-number="2717"></td>
        <td id="LC2717" class="blob-code js-file-line">    <span class="pl-c1">printImageInt</span>(meta,imageInt);</td>
      </tr>
      <tr>
        <td id="L2718" class="blob-num js-line-number" data-line-number="2718"></td>
        <td id="LC2718" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2719" class="blob-num js-line-number" data-line-number="2719"></td>
        <td id="LC2719" class="blob-code js-file-line">    <span class="pl-c1">free</span>(imageInt);</td>
      </tr>
      <tr>
        <td id="L2720" class="blob-num js-line-number" data-line-number="2720"></td>
        <td id="LC2720" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2721" class="blob-num js-line-number" data-line-number="2721"></td>
        <td id="LC2721" class="blob-code js-file-line">} <span class="pl-c">// printImageUChar</span></td>
      </tr>
      <tr>
        <td id="L2722" class="blob-num js-line-number" data-line-number="2722"></td>
        <td id="LC2722" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2723" class="blob-num js-line-number" data-line-number="2723"></td>
        <td id="LC2723" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2724" class="blob-num js-line-number" data-line-number="2724"></td>
        <td id="LC2724" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2725" class="blob-num js-line-number" data-line-number="2725"></td>
        <td id="LC2725" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2726" class="blob-num js-line-number" data-line-number="2726"></td>
        <td id="LC2726" class="blob-code js-file-line"><span class="pl-k">void</span> <span class="pl-en">vol2birdPrintIndexArrays</span>(<span class="pl-k">void</span>) {</td>
      </tr>
      <tr>
        <td id="L2727" class="blob-num js-line-number" data-line-number="2727"></td>
        <td id="LC2727" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2728" class="blob-num js-line-number" data-line-number="2728"></td>
        <td id="LC2728" class="blob-code js-file-line">    <span class="pl-k">if</span> (initializationSuccessful==<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L2729" class="blob-num js-line-number" data-line-number="2729"></td>
        <td id="LC2729" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>You need to initialize vol2bird before you can use it. Aborting.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2730" class="blob-num js-line-number" data-line-number="2730"></td>
        <td id="LC2730" class="blob-code js-file-line">        <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2731" class="blob-num js-line-number" data-line-number="2731"></td>
        <td id="LC2731" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2732" class="blob-num js-line-number" data-line-number="2732"></td>
        <td id="LC2732" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2733" class="blob-num js-line-number" data-line-number="2733"></td>
        <td id="LC2733" class="blob-code js-file-line">    <span class="pl-k">int</span> iLayer;</td>
      </tr>
      <tr>
        <td id="L2734" class="blob-num js-line-number" data-line-number="2734"></td>
        <td id="LC2734" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2735" class="blob-num js-line-number" data-line-number="2735"></td>
        <td id="LC2735" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>iLayer  iFrom   iTo     iTo-iFrom nWritten<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2736" class="blob-num js-line-number" data-line-number="2736"></td>
        <td id="LC2736" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L2737" class="blob-num js-line-number" data-line-number="2737"></td>
        <td id="LC2737" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%7d</span> <span class="pl-c1">%7d</span> <span class="pl-c1">%7d</span> <span class="pl-c1">%10d</span> <span class="pl-c1">%8d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L2738" class="blob-num js-line-number" data-line-number="2738"></td>
        <td id="LC2738" class="blob-code js-file-line">            iLayer, </td>
      </tr>
      <tr>
        <td id="L2739" class="blob-num js-line-number" data-line-number="2739"></td>
        <td id="LC2739" class="blob-code js-file-line">            indexFrom[iLayer], </td>
      </tr>
      <tr>
        <td id="L2740" class="blob-num js-line-number" data-line-number="2740"></td>
        <td id="LC2740" class="blob-code js-file-line">            indexTo[iLayer], </td>
      </tr>
      <tr>
        <td id="L2741" class="blob-num js-line-number" data-line-number="2741"></td>
        <td id="LC2741" class="blob-code js-file-line">            indexTo[iLayer] - indexFrom[iLayer], </td>
      </tr>
      <tr>
        <td id="L2742" class="blob-num js-line-number" data-line-number="2742"></td>
        <td id="LC2742" class="blob-code js-file-line">            nPointsWritten[iLayer]);</td>
      </tr>
      <tr>
        <td id="L2743" class="blob-num js-line-number" data-line-number="2743"></td>
        <td id="LC2743" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2744" class="blob-num js-line-number" data-line-number="2744"></td>
        <td id="LC2744" class="blob-code js-file-line">} <span class="pl-c">// vol2birdPrintIndexArrays</span></td>
      </tr>
      <tr>
        <td id="L2745" class="blob-num js-line-number" data-line-number="2745"></td>
        <td id="LC2745" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2746" class="blob-num js-line-number" data-line-number="2746"></td>
        <td id="LC2746" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2747" class="blob-num js-line-number" data-line-number="2747"></td>
        <td id="LC2747" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2748" class="blob-num js-line-number" data-line-number="2748"></td>
        <td id="LC2748" class="blob-code js-file-line"><span class="pl-k">void</span> <span class="pl-en">vol2birdPrintOptions</span>(<span class="pl-k">void</span>) {</td>
      </tr>
      <tr>
        <td id="L2749" class="blob-num js-line-number" data-line-number="2749"></td>
        <td id="LC2749" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2750" class="blob-num js-line-number" data-line-number="2750"></td>
        <td id="LC2750" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2751" class="blob-num js-line-number" data-line-number="2751"></td>
        <td id="LC2751" class="blob-code js-file-line">    <span class="pl-c">// this function prints vol2bird&#39;s configuration to stderr //</span></td>
      </tr>
      <tr>
        <td id="L2752" class="blob-num js-line-number" data-line-number="2752"></td>
        <td id="LC2752" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2753" class="blob-num js-line-number" data-line-number="2753"></td>
        <td id="LC2753" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2754" class="blob-num js-line-number" data-line-number="2754"></td>
        <td id="LC2754" class="blob-code js-file-line">    <span class="pl-k">if</span> (initializationSuccessful==<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L2755" class="blob-num js-line-number" data-line-number="2755"></td>
        <td id="LC2755" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>You need to initialize vol2bird before you can use it. Aborting.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2756" class="blob-num js-line-number" data-line-number="2756"></td>
        <td id="LC2756" class="blob-code js-file-line">        <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2757" class="blob-num js-line-number" data-line-number="2757"></td>
        <td id="LC2757" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2758" class="blob-num js-line-number" data-line-number="2758"></td>
        <td id="LC2758" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2759" class="blob-num js-line-number" data-line-number="2759"></td>
        <td id="LC2759" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n\n</span>vol2bird configuration:<span class="pl-cce">\n\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2760" class="blob-num js-line-number" data-line-number="2760"></td>
        <td id="LC2760" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2761" class="blob-num js-line-number" data-line-number="2761"></td>
        <td id="LC2761" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>absVDifMax<span class="pl-pds">&quot;</span></span>,absVDifMax);</td>
      </tr>
      <tr>
        <td id="L2762" class="blob-num js-line-number" data-line-number="2762"></td>
        <td id="LC2762" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>azimMax<span class="pl-pds">&quot;</span></span>,azimMax);</td>
      </tr>
      <tr>
        <td id="L2763" class="blob-num js-line-number" data-line-number="2763"></td>
        <td id="LC2763" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>azimMin<span class="pl-pds">&quot;</span></span>,azimMin);</td>
      </tr>
      <tr>
        <td id="L2764" class="blob-num js-line-number" data-line-number="2764"></td>
        <td id="LC2764" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>birdRadarCrossSection<span class="pl-pds">&quot;</span></span>,birdRadarCrossSection);</td>
      </tr>
      <tr>
        <td id="L2765" class="blob-num js-line-number" data-line-number="2765"></td>
        <td id="LC2765" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>cellClutterFractionMax<span class="pl-pds">&quot;</span></span>,cellClutterFractionMax);</td>
      </tr>
      <tr>
        <td id="L2766" class="blob-num js-line-number" data-line-number="2766"></td>
        <td id="LC2766" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>cellDbzMin<span class="pl-pds">&quot;</span></span>,cellDbzMin);</td>
      </tr>
      <tr>
        <td id="L2767" class="blob-num js-line-number" data-line-number="2767"></td>
        <td id="LC2767" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>cellStdDevMax<span class="pl-pds">&quot;</span></span>,cellStdDevMax);</td>
      </tr>
      <tr>
        <td id="L2768" class="blob-num js-line-number" data-line-number="2768"></td>
        <td id="LC2768" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>chisqMin<span class="pl-pds">&quot;</span></span>,chisqMin);</td>
      </tr>
      <tr>
        <td id="L2769" class="blob-num js-line-number" data-line-number="2769"></td>
        <td id="LC2769" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>clutterValueMin<span class="pl-pds">&quot;</span></span>,clutterValueMin);</td>
      </tr>
      <tr>
        <td id="L2770" class="blob-num js-line-number" data-line-number="2770"></td>
        <td id="LC2770" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>dbzFactor<span class="pl-pds">&quot;</span></span>,dbzFactor);</td>
      </tr>
      <tr>
        <td id="L2771" class="blob-num js-line-number" data-line-number="2771"></td>
        <td id="LC2771" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>dbzMax<span class="pl-pds">&quot;</span></span>,dbzMax);</td>
      </tr>
      <tr>
        <td id="L2772" class="blob-num js-line-number" data-line-number="2772"></td>
        <td id="LC2772" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>dbzThresMin<span class="pl-pds">&quot;</span></span>,dbzThresMin);</td>
      </tr>
      <tr>
        <td id="L2773" class="blob-num js-line-number" data-line-number="2773"></td>
        <td id="LC2773" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>fringeDist<span class="pl-pds">&quot;</span></span>,fringeDist);</td>
      </tr>
      <tr>
        <td id="L2774" class="blob-num js-line-number" data-line-number="2774"></td>
        <td id="LC2774" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>layerThickness<span class="pl-pds">&quot;</span></span>,layerThickness);</td>
      </tr>
      <tr>
        <td id="L2775" class="blob-num js-line-number" data-line-number="2775"></td>
        <td id="LC2775" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>nGatesCellMin<span class="pl-pds">&quot;</span></span>,nGatesCellMin);</td>
      </tr>
      <tr>
        <td id="L2776" class="blob-num js-line-number" data-line-number="2776"></td>
        <td id="LC2776" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>nAzimNeighborhood<span class="pl-pds">&quot;</span></span>,nAzimNeighborhood);</td>
      </tr>
      <tr>
        <td id="L2777" class="blob-num js-line-number" data-line-number="2777"></td>
        <td id="LC2777" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>nBinsGap<span class="pl-pds">&quot;</span></span>,nBinsGap);</td>
      </tr>
      <tr>
        <td id="L2778" class="blob-num js-line-number" data-line-number="2778"></td>
        <td id="LC2778" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>nCountMin<span class="pl-pds">&quot;</span></span>,nCountMin);</td>
      </tr>
      <tr>
        <td id="L2779" class="blob-num js-line-number" data-line-number="2779"></td>
        <td id="LC2779" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>nLayers<span class="pl-pds">&quot;</span></span>,nLayers);</td>
      </tr>
      <tr>
        <td id="L2780" class="blob-num js-line-number" data-line-number="2780"></td>
        <td id="LC2780" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>nObsGapMin<span class="pl-pds">&quot;</span></span>,nObsGapMin);</td>
      </tr>
      <tr>
        <td id="L2781" class="blob-num js-line-number" data-line-number="2781"></td>
        <td id="LC2781" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>nPointsIncludedMin<span class="pl-pds">&quot;</span></span>,nPointsIncludedMin);</td>
      </tr>
      <tr>
        <td id="L2782" class="blob-num js-line-number" data-line-number="2782"></td>
        <td id="LC2782" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>nRangNeighborhood<span class="pl-pds">&quot;</span></span>,nRangNeighborhood);</td>
      </tr>
      <tr>
        <td id="L2783" class="blob-num js-line-number" data-line-number="2783"></td>
        <td id="LC2783" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>radarWavelength<span class="pl-pds">&quot;</span></span>,radarWavelength);</td>
      </tr>
      <tr>
        <td id="L2784" class="blob-num js-line-number" data-line-number="2784"></td>
        <td id="LC2784" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>rangeMax<span class="pl-pds">&quot;</span></span>,rangeMax);</td>
      </tr>
      <tr>
        <td id="L2785" class="blob-num js-line-number" data-line-number="2785"></td>
        <td id="LC2785" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>rangeMin<span class="pl-pds">&quot;</span></span>,rangeMin);</td>
      </tr>
      <tr>
        <td id="L2786" class="blob-num js-line-number" data-line-number="2786"></td>
        <td id="LC2786" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>rCellMax<span class="pl-pds">&quot;</span></span>,rCellMax);</td>
      </tr>
      <tr>
        <td id="L2787" class="blob-num js-line-number" data-line-number="2787"></td>
        <td id="LC2787" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>refracIndex<span class="pl-pds">&quot;</span></span>,refracIndex);</td>
      </tr>
      <tr>
        <td id="L2788" class="blob-num js-line-number" data-line-number="2788"></td>
        <td id="LC2788" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>stdDevMinBird<span class="pl-pds">&quot;</span></span>,stdDevMinBird);</td>
      </tr>
      <tr>
        <td id="L2789" class="blob-num js-line-number" data-line-number="2789"></td>
        <td id="LC2789" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%c</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>useStaticClutterData<span class="pl-pds">&quot;</span></span>,useStaticClutterData == <span class="pl-c1">TRUE</span> ? <span class="pl-s"><span class="pl-pds">&#39;</span>T<span class="pl-pds">&#39;</span></span> : <span class="pl-s"><span class="pl-pds">&#39;</span>F<span class="pl-pds">&#39;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2790" class="blob-num js-line-number" data-line-number="2790"></td>
        <td id="LC2790" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%-25s</span> = <span class="pl-c1">%f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>vradMin<span class="pl-pds">&quot;</span></span>,vradMin);</td>
      </tr>
      <tr>
        <td id="L2791" class="blob-num js-line-number" data-line-number="2791"></td>
        <td id="LC2791" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2792" class="blob-num js-line-number" data-line-number="2792"></td>
        <td id="LC2792" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2793" class="blob-num js-line-number" data-line-number="2793"></td>
        <td id="LC2793" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2794" class="blob-num js-line-number" data-line-number="2794"></td>
        <td id="LC2794" class="blob-code js-file-line">}  <span class="pl-c">// vol2birdPrintOptions</span></td>
      </tr>
      <tr>
        <td id="L2795" class="blob-num js-line-number" data-line-number="2795"></td>
        <td id="LC2795" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2796" class="blob-num js-line-number" data-line-number="2796"></td>
        <td id="LC2796" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2797" class="blob-num js-line-number" data-line-number="2797"></td>
        <td id="LC2797" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2798" class="blob-num js-line-number" data-line-number="2798"></td>
        <td id="LC2798" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2799" class="blob-num js-line-number" data-line-number="2799"></td>
        <td id="LC2799" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2800" class="blob-num js-line-number" data-line-number="2800"></td>
        <td id="LC2800" class="blob-code js-file-line"><span class="pl-k">void</span> <span class="pl-en">vol2birdPrintPointsArray</span>(<span class="pl-k">void</span>) {</td>
      </tr>
      <tr>
        <td id="L2801" class="blob-num js-line-number" data-line-number="2801"></td>
        <td id="LC2801" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2802" class="blob-num js-line-number" data-line-number="2802"></td>
        <td id="LC2802" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2803" class="blob-num js-line-number" data-line-number="2803"></td>
        <td id="LC2803" class="blob-code js-file-line">    <span class="pl-c">// this function prints the &#39;points&#39; array to stderr //</span></td>
      </tr>
      <tr>
        <td id="L2804" class="blob-num js-line-number" data-line-number="2804"></td>
        <td id="LC2804" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2805" class="blob-num js-line-number" data-line-number="2805"></td>
        <td id="LC2805" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2806" class="blob-num js-line-number" data-line-number="2806"></td>
        <td id="LC2806" class="blob-code js-file-line">    <span class="pl-k">if</span> (initializationSuccessful==<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L2807" class="blob-num js-line-number" data-line-number="2807"></td>
        <td id="LC2807" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>You need to initialize vol2bird before you can use it. Aborting.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2808" class="blob-num js-line-number" data-line-number="2808"></td>
        <td id="LC2808" class="blob-code js-file-line">        <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2809" class="blob-num js-line-number" data-line-number="2809"></td>
        <td id="LC2809" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2810" class="blob-num js-line-number" data-line-number="2810"></td>
        <td id="LC2810" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2811" class="blob-num js-line-number" data-line-number="2811"></td>
        <td id="LC2811" class="blob-code js-file-line">    <span class="pl-k">int</span> iPoint;</td>
      </tr>
      <tr>
        <td id="L2812" class="blob-num js-line-number" data-line-number="2812"></td>
        <td id="LC2812" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2813" class="blob-num js-line-number" data-line-number="2813"></td>
        <td id="LC2813" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>iPoint  azim    elev    dbz         vrad        cell    gateCode  flags     <span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2814" class="blob-num js-line-number" data-line-number="2814"></td>
        <td id="LC2814" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2815" class="blob-num js-line-number" data-line-number="2815"></td>
        <td id="LC2815" class="blob-code js-file-line">    <span class="pl-k">for</span> (iPoint = <span class="pl-c1">0</span>; iPoint &lt; nRowsPoints * nColsPoints; iPoint+=nColsPoints) {</td>
      </tr>
      <tr>
        <td id="L2816" class="blob-num js-line-number" data-line-number="2816"></td>
        <td id="LC2816" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2817" class="blob-num js-line-number" data-line-number="2817"></td>
        <td id="LC2817" class="blob-code js-file-line">            <span class="pl-k">char</span> gateCodeStr[<span class="pl-c1">10</span>];  <span class="pl-c">// 9 bits plus 1 position for the null character &#39;\0&#39;</span></td>
      </tr>
      <tr>
        <td id="L2818" class="blob-num js-line-number" data-line-number="2818"></td>
        <td id="LC2818" class="blob-code js-file-line">            </td>
      </tr>
      <tr>
        <td id="L2819" class="blob-num js-line-number" data-line-number="2819"></td>
        <td id="LC2819" class="blob-code js-file-line">            <span class="pl-c1">printGateCode</span>(&amp;gateCodeStr[<span class="pl-c1">0</span>], (<span class="pl-k">int</span>) points[iPoint + gateCodeCol]);</td>
      </tr>
      <tr>
        <td id="L2820" class="blob-num js-line-number" data-line-number="2820"></td>
        <td id="LC2820" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2821" class="blob-num js-line-number" data-line-number="2821"></td>
        <td id="LC2821" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>  <span class="pl-c1">%6d</span><span class="pl-pds">&quot;</span></span>,    iPoint/nColsPoints);</td>
      </tr>
      <tr>
        <td id="L2822" class="blob-num js-line-number" data-line-number="2822"></td>
        <td id="LC2822" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>  <span class="pl-c1">%6.2f</span><span class="pl-pds">&quot;</span></span>,  points[iPoint + azimAngleCol]);</td>
      </tr>
      <tr>
        <td id="L2823" class="blob-num js-line-number" data-line-number="2823"></td>
        <td id="LC2823" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>  <span class="pl-c1">%6.2f</span><span class="pl-pds">&quot;</span></span>,  points[iPoint + elevAngleCol]);</td>
      </tr>
      <tr>
        <td id="L2824" class="blob-num js-line-number" data-line-number="2824"></td>
        <td id="LC2824" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>  <span class="pl-c1">%10.2f</span><span class="pl-pds">&quot;</span></span>, points[iPoint + dbzValueCol]);</td>
      </tr>
      <tr>
        <td id="L2825" class="blob-num js-line-number" data-line-number="2825"></td>
        <td id="LC2825" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>  <span class="pl-c1">%10.2f</span><span class="pl-pds">&quot;</span></span>, points[iPoint + vradValueCol]);</td>
      </tr>
      <tr>
        <td id="L2826" class="blob-num js-line-number" data-line-number="2826"></td>
        <td id="LC2826" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>  <span class="pl-c1">%6.0f</span><span class="pl-pds">&quot;</span></span>,  points[iPoint + cellValueCol]);</td>
      </tr>
      <tr>
        <td id="L2827" class="blob-num js-line-number" data-line-number="2827"></td>
        <td id="LC2827" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>  <span class="pl-c1">%8.0f</span><span class="pl-pds">&quot;</span></span>,  points[iPoint + gateCodeCol]);</td>
      </tr>
      <tr>
        <td id="L2828" class="blob-num js-line-number" data-line-number="2828"></td>
        <td id="LC2828" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>  <span class="pl-c1">%12s</span><span class="pl-pds">&quot;</span></span>,   gateCodeStr);</td>
      </tr>
      <tr>
        <td id="L2829" class="blob-num js-line-number" data-line-number="2829"></td>
        <td id="LC2829" class="blob-code js-file-line">            <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2830" class="blob-num js-line-number" data-line-number="2830"></td>
        <td id="LC2830" class="blob-code js-file-line">    }    </td>
      </tr>
      <tr>
        <td id="L2831" class="blob-num js-line-number" data-line-number="2831"></td>
        <td id="LC2831" class="blob-code js-file-line">} <span class="pl-c">// vol2birdPrintPointsArray</span></td>
      </tr>
      <tr>
        <td id="L2832" class="blob-num js-line-number" data-line-number="2832"></td>
        <td id="LC2832" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2833" class="blob-num js-line-number" data-line-number="2833"></td>
        <td id="LC2833" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2834" class="blob-num js-line-number" data-line-number="2834"></td>
        <td id="LC2834" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2835" class="blob-num js-line-number" data-line-number="2835"></td>
        <td id="LC2835" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2836" class="blob-num js-line-number" data-line-number="2836"></td>
        <td id="LC2836" class="blob-code js-file-line"><span class="pl-k">void</span> <span class="pl-en">printProfile</span>(<span class="pl-k">void</span>) {</td>
      </tr>
      <tr>
        <td id="L2837" class="blob-num js-line-number" data-line-number="2837"></td>
        <td id="LC2837" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2838" class="blob-num js-line-number" data-line-number="2838"></td>
        <td id="LC2838" class="blob-code js-file-line">    <span class="pl-k">if</span> (initializationSuccessful==<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L2839" class="blob-num js-line-number" data-line-number="2839"></td>
        <td id="LC2839" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>You need to initialize vol2bird before you can use it. Aborting.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2840" class="blob-num js-line-number" data-line-number="2840"></td>
        <td id="LC2840" class="blob-code js-file-line">        <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L2841" class="blob-num js-line-number" data-line-number="2841"></td>
        <td id="LC2841" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2842" class="blob-num js-line-number" data-line-number="2842"></td>
        <td id="LC2842" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2843" class="blob-num js-line-number" data-line-number="2843"></td>
        <td id="LC2843" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-cce">\n\n</span>Profile type: <span class="pl-c1">%d</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,iProfileTypeLast);</td>
      </tr>
      <tr>
        <td id="L2844" class="blob-num js-line-number" data-line-number="2844"></td>
        <td id="LC2844" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2845" class="blob-num js-line-number" data-line-number="2845"></td>
        <td id="LC2845" class="blob-code js-file-line">    <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>altmin-altmax: [u         ,v         ,w         ]; <span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L2846" class="blob-num js-line-number" data-line-number="2846"></td>
        <td id="LC2846" class="blob-code js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>hSpeed  , hDir    , chi     , hasGap  , dbzAvg  ,<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L2847" class="blob-num js-line-number" data-line-number="2847"></td>
        <td id="LC2847" class="blob-code js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span> nPoints, eta         , rhobird <span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2848" class="blob-num js-line-number" data-line-number="2848"></td>
        <td id="LC2848" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2849" class="blob-num js-line-number" data-line-number="2849"></td>
        <td id="LC2849" class="blob-code js-file-line">    <span class="pl-k">int</span> iLayer;</td>
      </tr>
      <tr>
        <td id="L2850" class="blob-num js-line-number" data-line-number="2850"></td>
        <td id="LC2850" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2851" class="blob-num js-line-number" data-line-number="2851"></td>
        <td id="LC2851" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = nLayers - <span class="pl-c1">1</span>; iLayer &gt;= <span class="pl-c1">0</span>; iLayer--) {</td>
      </tr>
      <tr>
        <td id="L2852" class="blob-num js-line-number" data-line-number="2852"></td>
        <td id="LC2852" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L2853" class="blob-num js-line-number" data-line-number="2853"></td>
        <td id="LC2853" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%6.0f</span>-<span class="pl-c1">%-6.0f</span>: [<span class="pl-c1">%10.2f</span>,<span class="pl-c1">%10.2f</span>,<span class="pl-c1">%10.2f</span>]; <span class="pl-c1">%8.2f</span>, <span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L2854" class="blob-num js-line-number" data-line-number="2854"></td>
        <td id="LC2854" class="blob-code js-file-line">        <span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">%8.1f</span>, <span class="pl-c1">%8.1f</span>, <span class="pl-c1">%8c</span>, <span class="pl-c1">%8.2f</span>, <span class="pl-c1">%7.0f</span>, <span class="pl-c1">%12.2f</span>, <span class="pl-c1">%8.2f</span><span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L2855" class="blob-num js-line-number" data-line-number="2855"></td>
        <td id="LC2855" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">0</span>],</td>
      </tr>
      <tr>
        <td id="L2856" class="blob-num js-line-number" data-line-number="2856"></td>
        <td id="LC2856" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">1</span>],</td>
      </tr>
      <tr>
        <td id="L2857" class="blob-num js-line-number" data-line-number="2857"></td>
        <td id="LC2857" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">2</span>],</td>
      </tr>
      <tr>
        <td id="L2858" class="blob-num js-line-number" data-line-number="2858"></td>
        <td id="LC2858" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">3</span>],</td>
      </tr>
      <tr>
        <td id="L2859" class="blob-num js-line-number" data-line-number="2859"></td>
        <td id="LC2859" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">4</span>],</td>
      </tr>
      <tr>
        <td id="L2860" class="blob-num js-line-number" data-line-number="2860"></td>
        <td id="LC2860" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">5</span>],</td>
      </tr>
      <tr>
        <td id="L2861" class="blob-num js-line-number" data-line-number="2861"></td>
        <td id="LC2861" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">6</span>],</td>
      </tr>
      <tr>
        <td id="L2862" class="blob-num js-line-number" data-line-number="2862"></td>
        <td id="LC2862" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">7</span>],</td>
      </tr>
      <tr>
        <td id="L2863" class="blob-num js-line-number" data-line-number="2863"></td>
        <td id="LC2863" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">8</span>] == <span class="pl-c1">TRUE</span> ? <span class="pl-s"><span class="pl-pds">&#39;</span>T<span class="pl-pds">&#39;</span></span> : <span class="pl-s"><span class="pl-pds">&#39;</span>F<span class="pl-pds">&#39;</span></span>,</td>
      </tr>
      <tr>
        <td id="L2864" class="blob-num js-line-number" data-line-number="2864"></td>
        <td id="LC2864" class="blob-code js-file-line">        profile[iLayer * nColsProfile +  <span class="pl-c1">9</span>],</td>
      </tr>
      <tr>
        <td id="L2865" class="blob-num js-line-number" data-line-number="2865"></td>
        <td id="LC2865" class="blob-code js-file-line">        profile[iLayer * nColsProfile + <span class="pl-c1">10</span>],</td>
      </tr>
      <tr>
        <td id="L2866" class="blob-num js-line-number" data-line-number="2866"></td>
        <td id="LC2866" class="blob-code js-file-line">        profile[iLayer * nColsProfile + <span class="pl-c1">11</span>],</td>
      </tr>
      <tr>
        <td id="L2867" class="blob-num js-line-number" data-line-number="2867"></td>
        <td id="LC2867" class="blob-code js-file-line">        profile[iLayer * nColsProfile + <span class="pl-c1">12</span>]);</td>
      </tr>
      <tr>
        <td id="L2868" class="blob-num js-line-number" data-line-number="2868"></td>
        <td id="LC2868" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2869" class="blob-num js-line-number" data-line-number="2869"></td>
        <td id="LC2869" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2870" class="blob-num js-line-number" data-line-number="2870"></td>
        <td id="LC2870" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2871" class="blob-num js-line-number" data-line-number="2871"></td>
        <td id="LC2871" class="blob-code js-file-line">} <span class="pl-c">// printProfile()</span></td>
      </tr>
      <tr>
        <td id="L2872" class="blob-num js-line-number" data-line-number="2872"></td>
        <td id="LC2872" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2873" class="blob-num js-line-number" data-line-number="2873"></td>
        <td id="LC2873" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2874" class="blob-num js-line-number" data-line-number="2874"></td>
        <td id="LC2874" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2875" class="blob-num js-line-number" data-line-number="2875"></td>
        <td id="LC2875" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2876" class="blob-num js-line-number" data-line-number="2876"></td>
        <td id="LC2876" class="blob-code js-file-line"><span class="pl-k">int</span> <span class="pl-en">vol2birdSetUp</span>(PolarVolume_t* volume) {</td>
      </tr>
      <tr>
        <td id="L2877" class="blob-num js-line-number" data-line-number="2877"></td>
        <td id="LC2877" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2878" class="blob-num js-line-number" data-line-number="2878"></td>
        <td id="LC2878" class="blob-code js-file-line">    <span class="pl-k">if</span> (<span class="pl-c1">readUserConfigOptions</span>() != <span class="pl-c1">0</span>) {</td>
      </tr>
      <tr>
        <td id="L2879" class="blob-num js-line-number" data-line-number="2879"></td>
        <td id="LC2879" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr, <span class="pl-s"><span class="pl-pds">&quot;</span>An error occurred while reading the user configuration file &#39;options.conf&#39;.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2880" class="blob-num js-line-number" data-line-number="2880"></td>
        <td id="LC2880" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>; </td>
      </tr>
      <tr>
        <td id="L2881" class="blob-num js-line-number" data-line-number="2881"></td>
        <td id="LC2881" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2882" class="blob-num js-line-number" data-line-number="2882"></td>
        <td id="LC2882" class="blob-code js-file-line">   </td>
      </tr>
      <tr>
        <td id="L2883" class="blob-num js-line-number" data-line-number="2883"></td>
        <td id="LC2883" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2884" class="blob-num js-line-number" data-line-number="2884"></td>
        <td id="LC2884" class="blob-code js-file-line">    <span class="pl-c">//              vol2bird options from options.conf               //</span></td>
      </tr>
      <tr>
        <td id="L2885" class="blob-num js-line-number" data-line-number="2885"></td>
        <td id="LC2885" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2886" class="blob-num js-line-number" data-line-number="2886"></td>
        <td id="LC2886" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2887" class="blob-num js-line-number" data-line-number="2887"></td>
        <td id="LC2887" class="blob-code js-file-line">    azimMax = <span class="pl-c1">cfg_getfloat</span>(cfg, <span class="pl-s"><span class="pl-pds">&quot;</span>AZIMMAX<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2888" class="blob-num js-line-number" data-line-number="2888"></td>
        <td id="LC2888" class="blob-code js-file-line">    azimMin = <span class="pl-c1">cfg_getfloat</span>(cfg, <span class="pl-s"><span class="pl-pds">&quot;</span>AZIMMIN<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2889" class="blob-num js-line-number" data-line-number="2889"></td>
        <td id="LC2889" class="blob-code js-file-line">    layerThickness = <span class="pl-c1">cfg_getfloat</span>(cfg, <span class="pl-s"><span class="pl-pds">&quot;</span>HLAYER<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2890" class="blob-num js-line-number" data-line-number="2890"></td>
        <td id="LC2890" class="blob-code js-file-line">    nLayers = <span class="pl-c1">cfg_getint</span>(cfg, <span class="pl-s"><span class="pl-pds">&quot;</span>NLAYER<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2891" class="blob-num js-line-number" data-line-number="2891"></td>
        <td id="LC2891" class="blob-code js-file-line">    rangeMax = <span class="pl-c1">cfg_getfloat</span>(cfg, <span class="pl-s"><span class="pl-pds">&quot;</span>RANGEMAX<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2892" class="blob-num js-line-number" data-line-number="2892"></td>
        <td id="LC2892" class="blob-code js-file-line">    rangeMin = <span class="pl-c1">cfg_getfloat</span>(cfg, <span class="pl-s"><span class="pl-pds">&quot;</span>RANGEMIN<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2893" class="blob-num js-line-number" data-line-number="2893"></td>
        <td id="LC2893" class="blob-code js-file-line">    radarWavelength = <span class="pl-c1">cfg_getfloat</span>(cfg, <span class="pl-s"><span class="pl-pds">&quot;</span>RADAR_WAVELENGTH_CM<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2894" class="blob-num js-line-number" data-line-number="2894"></td>
        <td id="LC2894" class="blob-code js-file-line">    useStaticClutterData = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>USE_STATIC_CLUTTER_DATA<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2895" class="blob-num js-line-number" data-line-number="2895"></td>
        <td id="LC2895" class="blob-code js-file-line">    printDbz = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_DBZ<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2896" class="blob-num js-line-number" data-line-number="2896"></td>
        <td id="LC2896" class="blob-code js-file-line">    printVrad = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_VRAD<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2897" class="blob-num js-line-number" data-line-number="2897"></td>
        <td id="LC2897" class="blob-code js-file-line">    printTex = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_TEXTURE<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2898" class="blob-num js-line-number" data-line-number="2898"></td>
        <td id="LC2898" class="blob-code js-file-line">    printCell = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_CELL<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2899" class="blob-num js-line-number" data-line-number="2899"></td>
        <td id="LC2899" class="blob-code js-file-line">    printCellProp = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_CELL_PROP<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2900" class="blob-num js-line-number" data-line-number="2900"></td>
        <td id="LC2900" class="blob-code js-file-line">    printClut = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_CLUT<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2901" class="blob-num js-line-number" data-line-number="2901"></td>
        <td id="LC2901" class="blob-code js-file-line">    printOptions = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_OPTIONS<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2902" class="blob-num js-line-number" data-line-number="2902"></td>
        <td id="LC2902" class="blob-code js-file-line">    printProfileVar = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_PROFILE<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2903" class="blob-num js-line-number" data-line-number="2903"></td>
        <td id="LC2903" class="blob-code js-file-line">    printPointsArray = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>PRINT_POINTS_ARRAY<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2904" class="blob-num js-line-number" data-line-number="2904"></td>
        <td id="LC2904" class="blob-code js-file-line">    fitVrad = <span class="pl-c1">cfg_getbool</span>(cfg,<span class="pl-s"><span class="pl-pds">&quot;</span>FIT_VRAD<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2905" class="blob-num js-line-number" data-line-number="2905"></td>
        <td id="LC2905" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2906" class="blob-num js-line-number" data-line-number="2906"></td>
        <td id="LC2906" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2907" class="blob-num js-line-number" data-line-number="2907"></td>
        <td id="LC2907" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2908" class="blob-num js-line-number" data-line-number="2908"></td>
        <td id="LC2908" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2909" class="blob-num js-line-number" data-line-number="2909"></td>
        <td id="LC2909" class="blob-code js-file-line">    <span class="pl-c">//              vol2bird options from constants.h                //</span></td>
      </tr>
      <tr>
        <td id="L2910" class="blob-num js-line-number" data-line-number="2910"></td>
        <td id="LC2910" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2911" class="blob-num js-line-number" data-line-number="2911"></td>
        <td id="LC2911" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2912" class="blob-num js-line-number" data-line-number="2912"></td>
        <td id="LC2912" class="blob-code js-file-line">    nGatesCellMin = AREACELL;</td>
      </tr>
      <tr>
        <td id="L2913" class="blob-num js-line-number" data-line-number="2913"></td>
        <td id="LC2913" class="blob-code js-file-line">    cellClutterFractionMax = CLUTPERCCELL;</td>
      </tr>
      <tr>
        <td id="L2914" class="blob-num js-line-number" data-line-number="2914"></td>
        <td id="LC2914" class="blob-code js-file-line">    cellDbzMin = DBZCELL;</td>
      </tr>
      <tr>
        <td id="L2915" class="blob-num js-line-number" data-line-number="2915"></td>
        <td id="LC2915" class="blob-code js-file-line">    chisqMin = CHISQMIN;</td>
      </tr>
      <tr>
        <td id="L2916" class="blob-num js-line-number" data-line-number="2916"></td>
        <td id="LC2916" class="blob-code js-file-line">    clutterValueMin = DBZCLUTTER;</td>
      </tr>
      <tr>
        <td id="L2917" class="blob-num js-line-number" data-line-number="2917"></td>
        <td id="LC2917" class="blob-code js-file-line">    dbzMax = DBZMAX;</td>
      </tr>
      <tr>
        <td id="L2918" class="blob-num js-line-number" data-line-number="2918"></td>
        <td id="LC2918" class="blob-code js-file-line">    dbzThresMin = DBZMIN;</td>
      </tr>
      <tr>
        <td id="L2919" class="blob-num js-line-number" data-line-number="2919"></td>
        <td id="LC2919" class="blob-code js-file-line">    fringeDist = FRINGEDIST;</td>
      </tr>
      <tr>
        <td id="L2920" class="blob-num js-line-number" data-line-number="2920"></td>
        <td id="LC2920" class="blob-code js-file-line">    nBinsGap = NBINSGAP;</td>
      </tr>
      <tr>
        <td id="L2921" class="blob-num js-line-number" data-line-number="2921"></td>
        <td id="LC2921" class="blob-code js-file-line">    nPointsIncludedMin = NDBZMIN;</td>
      </tr>
      <tr>
        <td id="L2922" class="blob-num js-line-number" data-line-number="2922"></td>
        <td id="LC2922" class="blob-code js-file-line">    nNeighborsMin = NEIGHBORS;</td>
      </tr>
      <tr>
        <td id="L2923" class="blob-num js-line-number" data-line-number="2923"></td>
        <td id="LC2923" class="blob-code js-file-line">    nObsGapMin = NOBSGAPMIN;</td>
      </tr>
      <tr>
        <td id="L2924" class="blob-num js-line-number" data-line-number="2924"></td>
        <td id="LC2924" class="blob-code js-file-line">    nAzimNeighborhood = NTEXBINAZIM;</td>
      </tr>
      <tr>
        <td id="L2925" class="blob-num js-line-number" data-line-number="2925"></td>
        <td id="LC2925" class="blob-code js-file-line">    nRangNeighborhood = NTEXBINRANG;</td>
      </tr>
      <tr>
        <td id="L2926" class="blob-num js-line-number" data-line-number="2926"></td>
        <td id="LC2926" class="blob-code js-file-line">    nCountMin = NTEXMIN; </td>
      </tr>
      <tr>
        <td id="L2927" class="blob-num js-line-number" data-line-number="2927"></td>
        <td id="LC2927" class="blob-code js-file-line">    refracIndex = REFRACTIVE_INDEX_OF_WATER;</td>
      </tr>
      <tr>
        <td id="L2928" class="blob-num js-line-number" data-line-number="2928"></td>
        <td id="LC2928" class="blob-code js-file-line">    birdRadarCrossSection = SIGMABIRD;</td>
      </tr>
      <tr>
        <td id="L2929" class="blob-num js-line-number" data-line-number="2929"></td>
        <td id="LC2929" class="blob-code js-file-line">    cellStdDevMax = STDEVCELL;</td>
      </tr>
      <tr>
        <td id="L2930" class="blob-num js-line-number" data-line-number="2930"></td>
        <td id="LC2930" class="blob-code js-file-line">    stdDevMinBird = STDEVBIRD;</td>
      </tr>
      <tr>
        <td id="L2931" class="blob-num js-line-number" data-line-number="2931"></td>
        <td id="LC2931" class="blob-code js-file-line">    absVDifMax = VDIFMAX;</td>
      </tr>
      <tr>
        <td id="L2932" class="blob-num js-line-number" data-line-number="2932"></td>
        <td id="LC2932" class="blob-code js-file-line">    vradMin = VRADMIN;</td>
      </tr>
      <tr>
        <td id="L2933" class="blob-num js-line-number" data-line-number="2933"></td>
        <td id="LC2933" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2934" class="blob-num js-line-number" data-line-number="2934"></td>
        <td id="LC2934" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2935" class="blob-num js-line-number" data-line-number="2935"></td>
        <td id="LC2935" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2936" class="blob-num js-line-number" data-line-number="2936"></td>
        <td id="LC2936" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2937" class="blob-num js-line-number" data-line-number="2937"></td>
        <td id="LC2937" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2938" class="blob-num js-line-number" data-line-number="2938"></td>
        <td id="LC2938" class="blob-code js-file-line">    <span class="pl-c">//                       some other variables                    //</span></td>
      </tr>
      <tr>
        <td id="L2939" class="blob-num js-line-number" data-line-number="2939"></td>
        <td id="LC2939" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2940" class="blob-num js-line-number" data-line-number="2940"></td>
        <td id="LC2940" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2941" class="blob-num js-line-number" data-line-number="2941"></td>
        <td id="LC2941" class="blob-code js-file-line">    rCellMax = rangeMax + <span class="pl-c1">5000</span>.<span class="pl-c1">0f</span>;</td>
      </tr>
      <tr>
        <td id="L2942" class="blob-num js-line-number" data-line-number="2942"></td>
        <td id="LC2942" class="blob-code js-file-line">    nDims = <span class="pl-c1">2</span>;</td>
      </tr>
      <tr>
        <td id="L2943" class="blob-num js-line-number" data-line-number="2943"></td>
        <td id="LC2943" class="blob-code js-file-line">    nParsFitted = <span class="pl-c1">3</span>;</td>
      </tr>
      <tr>
        <td id="L2944" class="blob-num js-line-number" data-line-number="2944"></td>
        <td id="LC2944" class="blob-code js-file-line">    dbzFactor = (<span class="pl-c1">pow</span>(refracIndex,<span class="pl-c1">2</span>) * <span class="pl-c1">1000</span> * <span class="pl-c1">pow</span>(PI,<span class="pl-c1">5</span>))/<span class="pl-c1">pow</span>(radarWavelength,<span class="pl-c1">4</span>);</td>
      </tr>
      <tr>
        <td id="L2945" class="blob-num js-line-number" data-line-number="2945"></td>
        <td id="LC2945" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2946" class="blob-num js-line-number" data-line-number="2946"></td>
        <td id="LC2946" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2947" class="blob-num js-line-number" data-line-number="2947"></td>
        <td id="LC2947" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2948" class="blob-num js-line-number" data-line-number="2948"></td>
        <td id="LC2948" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2949" class="blob-num js-line-number" data-line-number="2949"></td>
        <td id="LC2949" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2950" class="blob-num js-line-number" data-line-number="2950"></td>
        <td id="LC2950" class="blob-code js-file-line">    <span class="pl-c">//             lists of indices into the &#39;points&#39; array:         //</span></td>
      </tr>
      <tr>
        <td id="L2951" class="blob-num js-line-number" data-line-number="2951"></td>
        <td id="LC2951" class="blob-code js-file-line">    <span class="pl-c">//          where each altitude layer&#39;s data starts and ends     //</span></td>
      </tr>
      <tr>
        <td id="L2952" class="blob-num js-line-number" data-line-number="2952"></td>
        <td id="LC2952" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L2953" class="blob-num js-line-number" data-line-number="2953"></td>
        <td id="LC2953" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2954" class="blob-num js-line-number" data-line-number="2954"></td>
        <td id="LC2954" class="blob-code js-file-line">    <span class="pl-k">int</span> iLayer;</td>
      </tr>
      <tr>
        <td id="L2955" class="blob-num js-line-number" data-line-number="2955"></td>
        <td id="LC2955" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L2956" class="blob-num js-line-number" data-line-number="2956"></td>
        <td id="LC2956" class="blob-code js-file-line">    <span class="pl-c">// pre-allocate the list with start-from indexes for each </span></td>
      </tr>
      <tr>
        <td id="L2957" class="blob-num js-line-number" data-line-number="2957"></td>
        <td id="LC2957" class="blob-code js-file-line">    <span class="pl-c">// altitude bin in the profile</span></td>
      </tr>
      <tr>
        <td id="L2958" class="blob-num js-line-number" data-line-number="2958"></td>
        <td id="LC2958" class="blob-code js-file-line">    indexFrom = (<span class="pl-k">int</span>*) <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nLayers);</td>
      </tr>
      <tr>
        <td id="L2959" class="blob-num js-line-number" data-line-number="2959"></td>
        <td id="LC2959" class="blob-code js-file-line">    <span class="pl-k">if</span> (indexFrom == <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L2960" class="blob-num js-line-number" data-line-number="2960"></td>
        <td id="LC2960" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Error pre-allocating array &#39;indexFrom&#39;<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2961" class="blob-num js-line-number" data-line-number="2961"></td>
        <td id="LC2961" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2962" class="blob-num js-line-number" data-line-number="2962"></td>
        <td id="LC2962" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2963" class="blob-num js-line-number" data-line-number="2963"></td>
        <td id="LC2963" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L2964" class="blob-num js-line-number" data-line-number="2964"></td>
        <td id="LC2964" class="blob-code js-file-line">        indexFrom[iLayer] = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2965" class="blob-num js-line-number" data-line-number="2965"></td>
        <td id="LC2965" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2966" class="blob-num js-line-number" data-line-number="2966"></td>
        <td id="LC2966" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2967" class="blob-num js-line-number" data-line-number="2967"></td>
        <td id="LC2967" class="blob-code js-file-line">    <span class="pl-c">// pre-allocate the list with end-before indexes for each </span></td>
      </tr>
      <tr>
        <td id="L2968" class="blob-num js-line-number" data-line-number="2968"></td>
        <td id="LC2968" class="blob-code js-file-line">    <span class="pl-c">// altitude bin in the profile</span></td>
      </tr>
      <tr>
        <td id="L2969" class="blob-num js-line-number" data-line-number="2969"></td>
        <td id="LC2969" class="blob-code js-file-line">    indexTo = (<span class="pl-k">int</span>*) <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nLayers);</td>
      </tr>
      <tr>
        <td id="L2970" class="blob-num js-line-number" data-line-number="2970"></td>
        <td id="LC2970" class="blob-code js-file-line">    <span class="pl-k">if</span> (indexTo == <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L2971" class="blob-num js-line-number" data-line-number="2971"></td>
        <td id="LC2971" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Error pre-allocating array &#39;indexTo&#39;<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2972" class="blob-num js-line-number" data-line-number="2972"></td>
        <td id="LC2972" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2973" class="blob-num js-line-number" data-line-number="2973"></td>
        <td id="LC2973" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2974" class="blob-num js-line-number" data-line-number="2974"></td>
        <td id="LC2974" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L2975" class="blob-num js-line-number" data-line-number="2975"></td>
        <td id="LC2975" class="blob-code js-file-line">        indexTo[iLayer] = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L2976" class="blob-num js-line-number" data-line-number="2976"></td>
        <td id="LC2976" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2977" class="blob-num js-line-number" data-line-number="2977"></td>
        <td id="LC2977" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2978" class="blob-num js-line-number" data-line-number="2978"></td>
        <td id="LC2978" class="blob-code js-file-line">    <span class="pl-c">// pre-allocate the list containing TRUE or FALSE depending on the </span></td>
      </tr>
      <tr>
        <td id="L2979" class="blob-num js-line-number" data-line-number="2979"></td>
        <td id="LC2979" class="blob-code js-file-line">    <span class="pl-c">// results of calculating iProfileType == 3, which are needed when</span></td>
      </tr>
      <tr>
        <td id="L2980" class="blob-num js-line-number" data-line-number="2980"></td>
        <td id="LC2980" class="blob-code js-file-line">    <span class="pl-c">// calculating iProfileType == 1</span></td>
      </tr>
      <tr>
        <td id="L2981" class="blob-num js-line-number" data-line-number="2981"></td>
        <td id="LC2981" class="blob-code js-file-line">    scatterersAreNotBirds = (<span class="pl-k">int</span>*) <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nLayers);</td>
      </tr>
      <tr>
        <td id="L2982" class="blob-num js-line-number" data-line-number="2982"></td>
        <td id="LC2982" class="blob-code js-file-line">    <span class="pl-k">if</span> (scatterersAreNotBirds == <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L2983" class="blob-num js-line-number" data-line-number="2983"></td>
        <td id="LC2983" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Error pre-allocating array &#39;scatterersAreNotBirds&#39;<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2984" class="blob-num js-line-number" data-line-number="2984"></td>
        <td id="LC2984" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2985" class="blob-num js-line-number" data-line-number="2985"></td>
        <td id="LC2985" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2986" class="blob-num js-line-number" data-line-number="2986"></td>
        <td id="LC2986" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L2987" class="blob-num js-line-number" data-line-number="2987"></td>
        <td id="LC2987" class="blob-code js-file-line">        scatterersAreNotBirds[iLayer] = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2988" class="blob-num js-line-number" data-line-number="2988"></td>
        <td id="LC2988" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2989" class="blob-num js-line-number" data-line-number="2989"></td>
        <td id="LC2989" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2990" class="blob-num js-line-number" data-line-number="2990"></td>
        <td id="LC2990" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L2991" class="blob-num js-line-number" data-line-number="2991"></td>
        <td id="LC2991" class="blob-code js-file-line">    <span class="pl-c">// for each altitude layer, you need to remember how many points </span></td>
      </tr>
      <tr>
        <td id="L2992" class="blob-num js-line-number" data-line-number="2992"></td>
        <td id="LC2992" class="blob-code js-file-line">    <span class="pl-c">// were already written. This information is stored in the </span></td>
      </tr>
      <tr>
        <td id="L2993" class="blob-num js-line-number" data-line-number="2993"></td>
        <td id="LC2993" class="blob-code js-file-line">    <span class="pl-c">// &#39;nPointsWritten&#39; array</span></td>
      </tr>
      <tr>
        <td id="L2994" class="blob-num js-line-number" data-line-number="2994"></td>
        <td id="LC2994" class="blob-code js-file-line">    nPointsWritten = (<span class="pl-k">int</span>*) <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">int</span>) * nLayers);</td>
      </tr>
      <tr>
        <td id="L2995" class="blob-num js-line-number" data-line-number="2995"></td>
        <td id="LC2995" class="blob-code js-file-line">    <span class="pl-k">if</span> (nPointsWritten == <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L2996" class="blob-num js-line-number" data-line-number="2996"></td>
        <td id="LC2996" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Error pre-allocating array &#39;nPointsWritten&#39;<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L2997" class="blob-num js-line-number" data-line-number="2997"></td>
        <td id="LC2997" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L2998" class="blob-num js-line-number" data-line-number="2998"></td>
        <td id="LC2998" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L2999" class="blob-num js-line-number" data-line-number="2999"></td>
        <td id="LC2999" class="blob-code js-file-line">    <span class="pl-k">for</span> (iLayer = <span class="pl-c1">0</span>; iLayer &lt; nLayers; iLayer++) {</td>
      </tr>
      <tr>
        <td id="L3000" class="blob-num js-line-number" data-line-number="3000"></td>
        <td id="LC3000" class="blob-code js-file-line">        nPointsWritten[iLayer] = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L3001" class="blob-num js-line-number" data-line-number="3001"></td>
        <td id="LC3001" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L3002" class="blob-num js-line-number" data-line-number="3002"></td>
        <td id="LC3002" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3003" class="blob-num js-line-number" data-line-number="3003"></td>
        <td id="LC3003" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3004" class="blob-num js-line-number" data-line-number="3004"></td>
        <td id="LC3004" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3005" class="blob-num js-line-number" data-line-number="3005"></td>
        <td id="LC3005" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3006" class="blob-num js-line-number" data-line-number="3006"></td>
        <td id="LC3006" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L3007" class="blob-num js-line-number" data-line-number="3007"></td>
        <td id="LC3007" class="blob-code js-file-line">    <span class="pl-c">//               information about the &#39;points&#39; array            //</span></td>
      </tr>
      <tr>
        <td id="L3008" class="blob-num js-line-number" data-line-number="3008"></td>
        <td id="LC3008" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L3009" class="blob-num js-line-number" data-line-number="3009"></td>
        <td id="LC3009" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3010" class="blob-num js-line-number" data-line-number="3010"></td>
        <td id="LC3010" class="blob-code js-file-line">    nColsPoints = <span class="pl-c1">6</span>;</td>
      </tr>
      <tr>
        <td id="L3011" class="blob-num js-line-number" data-line-number="3011"></td>
        <td id="LC3011" class="blob-code js-file-line">    nRowsPoints = <span class="pl-c1">detSvdfitArraySize</span>(volume);</td>
      </tr>
      <tr>
        <td id="L3012" class="blob-num js-line-number" data-line-number="3012"></td>
        <td id="LC3012" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3013" class="blob-num js-line-number" data-line-number="3013"></td>
        <td id="LC3013" class="blob-code js-file-line">    azimAngleCol = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L3014" class="blob-num js-line-number" data-line-number="3014"></td>
        <td id="LC3014" class="blob-code js-file-line">    elevAngleCol = <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L3015" class="blob-num js-line-number" data-line-number="3015"></td>
        <td id="LC3015" class="blob-code js-file-line">    dbzValueCol = <span class="pl-c1">2</span>;</td>
      </tr>
      <tr>
        <td id="L3016" class="blob-num js-line-number" data-line-number="3016"></td>
        <td id="LC3016" class="blob-code js-file-line">    vradValueCol = <span class="pl-c1">3</span>;</td>
      </tr>
      <tr>
        <td id="L3017" class="blob-num js-line-number" data-line-number="3017"></td>
        <td id="LC3017" class="blob-code js-file-line">    cellValueCol = <span class="pl-c1">4</span>;</td>
      </tr>
      <tr>
        <td id="L3018" class="blob-num js-line-number" data-line-number="3018"></td>
        <td id="LC3018" class="blob-code js-file-line">    gateCodeCol = <span class="pl-c1">5</span>;</td>
      </tr>
      <tr>
        <td id="L3019" class="blob-num js-line-number" data-line-number="3019"></td>
        <td id="LC3019" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3020" class="blob-num js-line-number" data-line-number="3020"></td>
        <td id="LC3020" class="blob-code js-file-line">    <span class="pl-c">// pre-allocate the &#39;points&#39; array (note it has &#39;nColsPoints&#39;</span></td>
      </tr>
      <tr>
        <td id="L3021" class="blob-num js-line-number" data-line-number="3021"></td>
        <td id="LC3021" class="blob-code js-file-line">    <span class="pl-c">// pseudo-columns)</span></td>
      </tr>
      <tr>
        <td id="L3022" class="blob-num js-line-number" data-line-number="3022"></td>
        <td id="LC3022" class="blob-code js-file-line">    points = (<span class="pl-k">float</span>*) <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">float</span>) * nRowsPoints * nColsPoints);</td>
      </tr>
      <tr>
        <td id="L3023" class="blob-num js-line-number" data-line-number="3023"></td>
        <td id="LC3023" class="blob-code js-file-line">    <span class="pl-k">if</span> (points == <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L3024" class="blob-num js-line-number" data-line-number="3024"></td>
        <td id="LC3024" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Error pre-allocating array &#39;points&#39;.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>); </td>
      </tr>
      <tr>
        <td id="L3025" class="blob-num js-line-number" data-line-number="3025"></td>
        <td id="LC3025" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L3026" class="blob-num js-line-number" data-line-number="3026"></td>
        <td id="LC3026" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L3027" class="blob-num js-line-number" data-line-number="3027"></td>
        <td id="LC3027" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3028" class="blob-num js-line-number" data-line-number="3028"></td>
        <td id="LC3028" class="blob-code js-file-line">    <span class="pl-k">int</span> iRowPoints;</td>
      </tr>
      <tr>
        <td id="L3029" class="blob-num js-line-number" data-line-number="3029"></td>
        <td id="LC3029" class="blob-code js-file-line">    <span class="pl-k">int</span> iColPoints;</td>
      </tr>
      <tr>
        <td id="L3030" class="blob-num js-line-number" data-line-number="3030"></td>
        <td id="LC3030" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L3031" class="blob-num js-line-number" data-line-number="3031"></td>
        <td id="LC3031" class="blob-code js-file-line">    <span class="pl-k">for</span> (iRowPoints = <span class="pl-c1">0</span>; iRowPoints &lt; nRowsPoints; iRowPoints++) {</td>
      </tr>
      <tr>
        <td id="L3032" class="blob-num js-line-number" data-line-number="3032"></td>
        <td id="LC3032" class="blob-code js-file-line">        <span class="pl-k">for</span> (iColPoints = <span class="pl-c1">0</span>; iColPoints &lt; nColsPoints; iColPoints++) {</td>
      </tr>
      <tr>
        <td id="L3033" class="blob-num js-line-number" data-line-number="3033"></td>
        <td id="LC3033" class="blob-code js-file-line">            points[iRowPoints*nColsPoints + iColPoints] = NAN;</td>
      </tr>
      <tr>
        <td id="L3034" class="blob-num js-line-number" data-line-number="3034"></td>
        <td id="LC3034" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L3035" class="blob-num js-line-number" data-line-number="3035"></td>
        <td id="LC3035" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L3036" class="blob-num js-line-number" data-line-number="3036"></td>
        <td id="LC3036" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3037" class="blob-num js-line-number" data-line-number="3037"></td>
        <td id="LC3037" class="blob-code js-file-line">    <span class="pl-c">// information about the flagfields of &#39;gateCode&#39;</span></td>
      </tr>
      <tr>
        <td id="L3038" class="blob-num js-line-number" data-line-number="3038"></td>
        <td id="LC3038" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L3039" class="blob-num js-line-number" data-line-number="3039"></td>
        <td id="LC3039" class="blob-code js-file-line">    flagPositionStaticClutter = <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L3040" class="blob-num js-line-number" data-line-number="3040"></td>
        <td id="LC3040" class="blob-code js-file-line">    flagPositionDynamicClutter = <span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L3041" class="blob-num js-line-number" data-line-number="3041"></td>
        <td id="LC3041" class="blob-code js-file-line">    flagPositionDynamicClutterFringe = <span class="pl-c1">2</span>;</td>
      </tr>
      <tr>
        <td id="L3042" class="blob-num js-line-number" data-line-number="3042"></td>
        <td id="LC3042" class="blob-code js-file-line">    flagPositionVradMissing = <span class="pl-c1">3</span>;</td>
      </tr>
      <tr>
        <td id="L3043" class="blob-num js-line-number" data-line-number="3043"></td>
        <td id="LC3043" class="blob-code js-file-line">    flagPositionDbzTooHighForBirds = <span class="pl-c1">4</span>;</td>
      </tr>
      <tr>
        <td id="L3044" class="blob-num js-line-number" data-line-number="3044"></td>
        <td id="LC3044" class="blob-code js-file-line">    flagPositionVradTooLow = <span class="pl-c1">5</span>;</td>
      </tr>
      <tr>
        <td id="L3045" class="blob-num js-line-number" data-line-number="3045"></td>
        <td id="LC3045" class="blob-code js-file-line">    flagPositionVDifMax = <span class="pl-c1">6</span>;</td>
      </tr>
      <tr>
        <td id="L3046" class="blob-num js-line-number" data-line-number="3046"></td>
        <td id="LC3046" class="blob-code js-file-line">    flagPositionAzimTooLow = <span class="pl-c1">7</span>;</td>
      </tr>
      <tr>
        <td id="L3047" class="blob-num js-line-number" data-line-number="3047"></td>
        <td id="LC3047" class="blob-code js-file-line">    flagPositionAzimTooHigh = <span class="pl-c1">8</span>;</td>
      </tr>
      <tr>
        <td id="L3048" class="blob-num js-line-number" data-line-number="3048"></td>
        <td id="LC3048" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3049" class="blob-num js-line-number" data-line-number="3049"></td>
        <td id="LC3049" class="blob-code js-file-line">    <span class="pl-c">// construct the &#39;points&#39; array</span></td>
      </tr>
      <tr>
        <td id="L3050" class="blob-num js-line-number" data-line-number="3050"></td>
        <td id="LC3050" class="blob-code js-file-line">    <span class="pl-c1">constructPointsArray</span>(volume);</td>
      </tr>
      <tr>
        <td id="L3051" class="blob-num js-line-number" data-line-number="3051"></td>
        <td id="LC3051" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3052" class="blob-num js-line-number" data-line-number="3052"></td>
        <td id="LC3052" class="blob-code js-file-line">    <span class="pl-c">// classify the gates based on the data in &#39;points&#39;</span></td>
      </tr>
      <tr>
        <td id="L3053" class="blob-num js-line-number" data-line-number="3053"></td>
        <td id="LC3053" class="blob-code js-file-line">    <span class="pl-c1">classifyGatesSimple</span>();</td>
      </tr>
      <tr>
        <td id="L3054" class="blob-num js-line-number" data-line-number="3054"></td>
        <td id="LC3054" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3055" class="blob-num js-line-number" data-line-number="3055"></td>
        <td id="LC3055" class="blob-code js-file-line">    <span class="pl-k">if</span> (printPointsArray == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L3056" class="blob-num js-line-number" data-line-number="3056"></td>
        <td id="LC3056" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3057" class="blob-num js-line-number" data-line-number="3057"></td>
        <td id="LC3057" class="blob-code js-file-line">        <span class="pl-c1">vol2birdPrintIndexArrays</span>();</td>
      </tr>
      <tr>
        <td id="L3058" class="blob-num js-line-number" data-line-number="3058"></td>
        <td id="LC3058" class="blob-code js-file-line">        <span class="pl-c1">vol2birdPrintPointsArray</span>();</td>
      </tr>
      <tr>
        <td id="L3059" class="blob-num js-line-number" data-line-number="3059"></td>
        <td id="LC3059" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3060" class="blob-num js-line-number" data-line-number="3060"></td>
        <td id="LC3060" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L3061" class="blob-num js-line-number" data-line-number="3061"></td>
        <td id="LC3061" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3062" class="blob-num js-line-number" data-line-number="3062"></td>
        <td id="LC3062" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L3063" class="blob-num js-line-number" data-line-number="3063"></td>
        <td id="LC3063" class="blob-code js-file-line">    <span class="pl-c">//              information about the &#39;profile&#39; array            //</span></td>
      </tr>
      <tr>
        <td id="L3064" class="blob-num js-line-number" data-line-number="3064"></td>
        <td id="LC3064" class="blob-code js-file-line">    <span class="pl-c">// ------------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L3065" class="blob-num js-line-number" data-line-number="3065"></td>
        <td id="LC3065" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3066" class="blob-num js-line-number" data-line-number="3066"></td>
        <td id="LC3066" class="blob-code js-file-line">    nProfileTypes = <span class="pl-c1">3</span>;</td>
      </tr>
      <tr>
        <td id="L3067" class="blob-num js-line-number" data-line-number="3067"></td>
        <td id="LC3067" class="blob-code js-file-line">    nRowsProfile = nLayers;</td>
      </tr>
      <tr>
        <td id="L3068" class="blob-num js-line-number" data-line-number="3068"></td>
        <td id="LC3068" class="blob-code js-file-line">    nColsProfile = <span class="pl-c1">13</span>; </td>
      </tr>
      <tr>
        <td id="L3069" class="blob-num js-line-number" data-line-number="3069"></td>
        <td id="LC3069" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L3070" class="blob-num js-line-number" data-line-number="3070"></td>
        <td id="LC3070" class="blob-code js-file-line">    <span class="pl-c">// pre-allocate the array holding any profiled data (note it has </span></td>
      </tr>
      <tr>
        <td id="L3071" class="blob-num js-line-number" data-line-number="3071"></td>
        <td id="LC3071" class="blob-code js-file-line">    <span class="pl-c">// &#39;nColsProfile&#39; pseudocolumns):</span></td>
      </tr>
      <tr>
        <td id="L3072" class="blob-num js-line-number" data-line-number="3072"></td>
        <td id="LC3072" class="blob-code js-file-line">    profile = (<span class="pl-k">float</span>*) <span class="pl-c1">malloc</span>(<span class="pl-k">sizeof</span>(<span class="pl-k">float</span>) * nRowsProfile * nColsProfile);</td>
      </tr>
      <tr>
        <td id="L3073" class="blob-num js-line-number" data-line-number="3073"></td>
        <td id="LC3073" class="blob-code js-file-line">    <span class="pl-k">if</span> (profile == <span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L3074" class="blob-num js-line-number" data-line-number="3074"></td>
        <td id="LC3074" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>Error pre-allocating array &#39;profile&#39;.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>); </td>
      </tr>
      <tr>
        <td id="L3075" class="blob-num js-line-number" data-line-number="3075"></td>
        <td id="LC3075" class="blob-code js-file-line">        <span class="pl-k">return</span> -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L3076" class="blob-num js-line-number" data-line-number="3076"></td>
        <td id="LC3076" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L3077" class="blob-num js-line-number" data-line-number="3077"></td>
        <td id="LC3077" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3078" class="blob-num js-line-number" data-line-number="3078"></td>
        <td id="LC3078" class="blob-code js-file-line">    <span class="pl-k">int</span> iRowProfile;</td>
      </tr>
      <tr>
        <td id="L3079" class="blob-num js-line-number" data-line-number="3079"></td>
        <td id="LC3079" class="blob-code js-file-line">    <span class="pl-k">int</span> iColProfile;</td>
      </tr>
      <tr>
        <td id="L3080" class="blob-num js-line-number" data-line-number="3080"></td>
        <td id="LC3080" class="blob-code js-file-line">        </td>
      </tr>
      <tr>
        <td id="L3081" class="blob-num js-line-number" data-line-number="3081"></td>
        <td id="LC3081" class="blob-code js-file-line">    <span class="pl-k">for</span> (iRowProfile = <span class="pl-c1">0</span>; iRowProfile &lt; nRowsProfile; iRowProfile++) {</td>
      </tr>
      <tr>
        <td id="L3082" class="blob-num js-line-number" data-line-number="3082"></td>
        <td id="LC3082" class="blob-code js-file-line">        <span class="pl-k">for</span> (iColProfile = <span class="pl-c1">0</span>; iColProfile &lt; nColsProfile; iColProfile++) {</td>
      </tr>
      <tr>
        <td id="L3083" class="blob-num js-line-number" data-line-number="3083"></td>
        <td id="LC3083" class="blob-code js-file-line">            profile[iRowProfile*nColsProfile + iColProfile] = NAN;</td>
      </tr>
      <tr>
        <td id="L3084" class="blob-num js-line-number" data-line-number="3084"></td>
        <td id="LC3084" class="blob-code js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L3085" class="blob-num js-line-number" data-line-number="3085"></td>
        <td id="LC3085" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L3086" class="blob-num js-line-number" data-line-number="3086"></td>
        <td id="LC3086" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3087" class="blob-num js-line-number" data-line-number="3087"></td>
        <td id="LC3087" class="blob-code js-file-line">    iProfileTypeLast = -<span class="pl-c1">1</span>;</td>
      </tr>
      <tr>
        <td id="L3088" class="blob-num js-line-number" data-line-number="3088"></td>
        <td id="LC3088" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3089" class="blob-num js-line-number" data-line-number="3089"></td>
        <td id="LC3089" class="blob-code js-file-line">    initializationSuccessful = <span class="pl-c1">TRUE</span>;</td>
      </tr>
      <tr>
        <td id="L3090" class="blob-num js-line-number" data-line-number="3090"></td>
        <td id="LC3090" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3091" class="blob-num js-line-number" data-line-number="3091"></td>
        <td id="LC3091" class="blob-code js-file-line">    <span class="pl-k">if</span> (printOptions == <span class="pl-c1">TRUE</span>) {</td>
      </tr>
      <tr>
        <td id="L3092" class="blob-num js-line-number" data-line-number="3092"></td>
        <td id="LC3092" class="blob-code js-file-line">        <span class="pl-c1">vol2birdPrintOptions</span>();</td>
      </tr>
      <tr>
        <td id="L3093" class="blob-num js-line-number" data-line-number="3093"></td>
        <td id="LC3093" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L3094" class="blob-num js-line-number" data-line-number="3094"></td>
        <td id="LC3094" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3095" class="blob-num js-line-number" data-line-number="3095"></td>
        <td id="LC3095" class="blob-code js-file-line">    <span class="pl-k">return</span> <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L3096" class="blob-num js-line-number" data-line-number="3096"></td>
        <td id="LC3096" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3097" class="blob-num js-line-number" data-line-number="3097"></td>
        <td id="LC3097" class="blob-code js-file-line">} <span class="pl-c">// vol2birdSetUp</span></td>
      </tr>
      <tr>
        <td id="L3098" class="blob-num js-line-number" data-line-number="3098"></td>
        <td id="LC3098" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3099" class="blob-num js-line-number" data-line-number="3099"></td>
        <td id="LC3099" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3100" class="blob-num js-line-number" data-line-number="3100"></td>
        <td id="LC3100" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3101" class="blob-num js-line-number" data-line-number="3101"></td>
        <td id="LC3101" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3102" class="blob-num js-line-number" data-line-number="3102"></td>
        <td id="LC3102" class="blob-code js-file-line"><span class="pl-k">void</span> <span class="pl-en">vol2birdTearDown</span>() {</td>
      </tr>
      <tr>
        <td id="L3103" class="blob-num js-line-number" data-line-number="3103"></td>
        <td id="LC3103" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L3104" class="blob-num js-line-number" data-line-number="3104"></td>
        <td id="LC3104" class="blob-code js-file-line">    <span class="pl-c">// ---------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L3105" class="blob-num js-line-number" data-line-number="3105"></td>
        <td id="LC3105" class="blob-code js-file-line">    <span class="pl-c">// free the memory that was previously allocated for vol2bird //</span></td>
      </tr>
      <tr>
        <td id="L3106" class="blob-num js-line-number" data-line-number="3106"></td>
        <td id="LC3106" class="blob-code js-file-line">    <span class="pl-c">// ---------------------------------------------------------- //</span></td>
      </tr>
      <tr>
        <td id="L3107" class="blob-num js-line-number" data-line-number="3107"></td>
        <td id="LC3107" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3108" class="blob-num js-line-number" data-line-number="3108"></td>
        <td id="LC3108" class="blob-code js-file-line">    <span class="pl-k">if</span> (initializationSuccessful==<span class="pl-c1">FALSE</span>) {</td>
      </tr>
      <tr>
        <td id="L3109" class="blob-num js-line-number" data-line-number="3109"></td>
        <td id="LC3109" class="blob-code js-file-line">        <span class="pl-c1">fprintf</span>(stderr,<span class="pl-s"><span class="pl-pds">&quot;</span>You need to initialize vol2bird before you can use it. Aborting.<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L3110" class="blob-num js-line-number" data-line-number="3110"></td>
        <td id="LC3110" class="blob-code js-file-line">        <span class="pl-k">return</span>;</td>
      </tr>
      <tr>
        <td id="L3111" class="blob-num js-line-number" data-line-number="3111"></td>
        <td id="LC3111" class="blob-code js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L3112" class="blob-num js-line-number" data-line-number="3112"></td>
        <td id="LC3112" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3113" class="blob-num js-line-number" data-line-number="3113"></td>
        <td id="LC3113" class="blob-code js-file-line">    <span class="pl-c">// free the points array, the indexes into it, the counters, as well</span></td>
      </tr>
      <tr>
        <td id="L3114" class="blob-num js-line-number" data-line-number="3114"></td>
        <td id="LC3114" class="blob-code js-file-line">    <span class="pl-c">// as the profile data array</span></td>
      </tr>
      <tr>
        <td id="L3115" class="blob-num js-line-number" data-line-number="3115"></td>
        <td id="LC3115" class="blob-code js-file-line">    <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) points);</td>
      </tr>
      <tr>
        <td id="L3116" class="blob-num js-line-number" data-line-number="3116"></td>
        <td id="LC3116" class="blob-code js-file-line">    <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) profile);</td>
      </tr>
      <tr>
        <td id="L3117" class="blob-num js-line-number" data-line-number="3117"></td>
        <td id="LC3117" class="blob-code js-file-line">    <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) indexFrom);</td>
      </tr>
      <tr>
        <td id="L3118" class="blob-num js-line-number" data-line-number="3118"></td>
        <td id="LC3118" class="blob-code js-file-line">    <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) indexTo);</td>
      </tr>
      <tr>
        <td id="L3119" class="blob-num js-line-number" data-line-number="3119"></td>
        <td id="LC3119" class="blob-code js-file-line">    <span class="pl-c1">free</span>((<span class="pl-k">void</span>*) nPointsWritten);</td>
      </tr>
      <tr>
        <td id="L3120" class="blob-num js-line-number" data-line-number="3120"></td>
        <td id="LC3120" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L3121" class="blob-num js-line-number" data-line-number="3121"></td>
        <td id="LC3121" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L3122" class="blob-num js-line-number" data-line-number="3122"></td>
        <td id="LC3122" class="blob-code js-file-line">    <span class="pl-c">// free the memory that holds the user configurable options</span></td>
      </tr>
      <tr>
        <td id="L3123" class="blob-num js-line-number" data-line-number="3123"></td>
        <td id="LC3123" class="blob-code js-file-line">    <span class="pl-c1">cfg_free</span>(cfg);</td>
      </tr>
      <tr>
        <td id="L3124" class="blob-num js-line-number" data-line-number="3124"></td>
        <td id="LC3124" class="blob-code js-file-line">    </td>
      </tr>
      <tr>
        <td id="L3125" class="blob-num js-line-number" data-line-number="3125"></td>
        <td id="LC3125" class="blob-code js-file-line">    <span class="pl-c">// reset this variable to its initial value</span></td>
      </tr>
      <tr>
        <td id="L3126" class="blob-num js-line-number" data-line-number="3126"></td>
        <td id="LC3126" class="blob-code js-file-line">    initializationSuccessful = <span class="pl-c1">FALSE</span>;</td>
      </tr>
      <tr>
        <td id="L3127" class="blob-num js-line-number" data-line-number="3127"></td>
        <td id="LC3127" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3128" class="blob-num js-line-number" data-line-number="3128"></td>
        <td id="LC3128" class="blob-code js-file-line">} <span class="pl-c">// vol2birdTearDown</span></td>
      </tr>
      <tr>
        <td id="L3129" class="blob-num js-line-number" data-line-number="3129"></td>
        <td id="LC3129" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3130" class="blob-num js-line-number" data-line-number="3130"></td>
        <td id="LC3130" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3131" class="blob-num js-line-number" data-line-number="3131"></td>
        <td id="LC3131" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3132" class="blob-num js-line-number" data-line-number="3132"></td>
        <td id="LC3132" class="blob-code js-file-line">
</td>
      </tr>
</table>

  </div>

</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" action="" class="js-jump-to-line-form" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="btn">Go</button>
</form></div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer" role="contentinfo">
    <ul class="site-footer-links right">
        <li><a href="https://status.github.com/" data-ga-click="Footer, go to status, text:status">Status</a></li>
      <li><a href="https://developer.github.com" data-ga-click="Footer, go to api, text:api">API</a></li>
      <li><a href="https://training.github.com" data-ga-click="Footer, go to training, text:training">Training</a></li>
      <li><a href="https://shop.github.com" data-ga-click="Footer, go to shop, text:shop">Shop</a></li>
        <li><a href="https://github.com/blog" data-ga-click="Footer, go to blog, text:blog">Blog</a></li>
        <li><a href="https://github.com/about" data-ga-click="Footer, go to about, text:about">About</a></li>

    </ul>

    <a href="https://github.com" aria-label="Homepage">
      <span class="mega-octicon octicon-mark-github" title="GitHub"></span>
</a>
    <ul class="site-footer-links">
      <li>&copy; 2015 <span title="0.10284s from github-fe129-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="https://github.com/site/terms" data-ga-click="Footer, go to terms, text:terms">Terms</a></li>
        <li><a href="https://github.com/site/privacy" data-ga-click="Footer, go to privacy, text:privacy">Privacy</a></li>
        <li><a href="https://github.com/security" data-ga-click="Footer, go to security, text:security">Security</a></li>
        <li><a href="https://github.com/contact" data-ga-click="Footer, go to contact, text:contact">Contact</a></li>
    </ul>
  </div>
</div>


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-suggester-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="fullscreen-contents js-fullscreen-contents" placeholder=""></textarea>
      <div class="suggester-container">
        <div class="suggester fullscreen-suggester js-suggester js-navigation-container"></div>
      </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped tooltipped-w" aria-label="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped tooltipped-w"
      aria-label="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    
    

    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-x flash-close js-ajax-error-dismiss" aria-label="Dismiss error"></a>
      Something went wrong with that request. Please try again.
    </div>


      <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/frameworks-d22b59d0085e83b7549ba4341ec9e68f80c2f29c8e49213ee182003dc8d568c6.js"></script>
      <script async="async" crossorigin="anonymous" src="https://assets-cdn.github.com/assets/github-3b4fde20974117506181132ac0ad3588b632eb7b8b8347002c439b67027b0227.js"></script>
      
      

  </body>
</html>


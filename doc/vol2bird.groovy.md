# Baltrad example script

To run vol2bird from Baltrad, you need to add a Route that runs vol2bird through RAVE.
The below script will run vol2bird if configured as a Script Route in the Baltrad-DEX
user interface. You'll need an adapter, and to schedule the route, as follows.

## Create an Adaptor

Go to Adaptors and click Create

Name: RAVE
Type: XMLRPC
URI: http://eslt0045:8085/RAVE
Timeout: 5000

## Create a Route

Go to Routes -> Create script

Name: vol2bird
Author: Jurriaan and Lourens
Active: yes
Description: Run vol2bird on incoming files
Recipients: RAVE (select)
Script: See vol2bird.groovy below

## Schedule the route

Go to Schedule and click Create
Schedule vol2bird every minute

## vol2bird.groovy

```
package eu.baltrad.vol2bird;

import java.util.HashSet;
import java.util.Set;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import eu.baltrad.bdb.db.FileEntry;
import eu.baltrad.bdb.oh5.Metadata;
import eu.baltrad.bdb.oh5.TemplateMetadataNamer;
import eu.baltrad.bdb.util.Date;
import eu.baltrad.bdb.util.Time;
import eu.baltrad.beast.ManagerContext;
import eu.baltrad.beast.db.Catalog;
import eu.baltrad.beast.db.CatalogEntry;
import eu.baltrad.beast.message.IBltMessage;
import eu.baltrad.beast.rules.IScriptableRule;

import eu.baltrad.beast.message.mo.BltDataMessage;
import eu.baltrad.beast.message.mo.BltGenerateMessage;

class Vol2Bird implements IScriptableRule {
  private static Set<String> SUPPORTED_RADARS = new HashSet<String>([
    "bejab", "bewid", "bezav", "dkbor", "dkrom",
    "dksin", "dkste", "dkvir", "eehar", "eesur",
    "fianj", "fiika", "fikes", "fikor", "fikuo",
    "filuo", "fiuta", "fivan", "fivim", "frabb",
    "frale", "frbla", "frbol", "frbor", "frbou",
    "frcae", "frche", "frcol", "frgre", "frlep",
    "frmcl", "frmom", "frmtc", "frnan", "frnim",
    "fropo", "frpla", "frtou", "frtra", "frtre",
    "hrbil", "hrosi", "nldbl", "nldhl", "searl",
    "sease", "sehud", "sekir", "sekkr", "selek",
    "selul", "seosu", "seovi", "sevar", "sevil",
    "silis", "sipas"
  ]);

  private Logger logger = LogManager.getLogger(Vol2Bird.class);

  public Vol2Bird() {
  }

  private boolean isSupported(String radar) {
    return SUPPORTED_RADARS.contains(radar);
  }

  @Override
  public IBltMessage handle(IBltMessage bltmsg) {
    BltGenerateMessage result = null;

    if (bltmsg != null && bltmsg instanceof BltDataMessage) {
      FileEntry fe = ((BltDataMessage)bltmsg).getFileEntry();

      String object = fe.getMetadata().getWhatObject();
      if (object != null && object.equals("PVOL")) {
        Metadata metadata = fe.getMetadata();
        String radar = metadata.getAttribute("/_bdb/source_name");

        if (isSupported(radar)) {
          Catalog cat = ManagerContext.getCatalog();

          result = new BltGenerateMessage();
          result.setAlgorithm("eu.baltrad.beast.vol2bird");
          result.setFiles([fe.getUuid().toString()] as String[])
          result.setArguments(new String[0]);
        }
        else {
          logger.error(sprintf("Radar %s not supported by vol2bird", radar));
          logger.info("Supported radars:");
          for (String s : SUPPORTED_RADARS) {
            logger.info(s);
          }
        }
      }
    }
    return result;
  }
}
```


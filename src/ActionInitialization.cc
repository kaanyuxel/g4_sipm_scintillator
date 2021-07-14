#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction());
  RunAction* theRunAction = new RunAction();
  SetUserAction(theRunAction);
  EventAction* eventAction = new EventAction();
 	SetUserAction(eventAction);
}

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction());
}

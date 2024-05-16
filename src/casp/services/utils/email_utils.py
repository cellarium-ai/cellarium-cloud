import logging
import typing as t

from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

logger = logging.getLogger(__name__)

TEMPLATE_DIR = "src/casp/services/admin/templates"


class EmailSender:

    def __init__(
        self,
        sendgrid_client: t.Optional[t.Union[SendGridAPIClient, str]],
        from_address: str,
        template_dir: str = TEMPLATE_DIR,
    ) -> None:
        """
        Initialize EmailSender class.

        :param sendgrid_client: SendGrid API key or client object. If this value is empty, the class will not send emails.
        :param from_address: Email address to send emails from.
        """
        if sendgrid_client is None:
            logger.warning("No SendGrid key provided.  Emails will not be sent.")
            self.sg = None
            return

        if isinstance(sendgrid_client, str):
            self.sg = SendGridAPIClient(sendgrid_client)
        elif isinstance(sendgrid_client, SendGridAPIClient):
            self.sg = sendgrid_client
        else:
            raise ValueError("sendgrid_client must be a SendGridAPIClient object or a string.")

        self.from_address = from_address
        self.template_dir = template_dir

    def _send_email(self, to: str, subject: str, content_path_html: str, content_path_plain: str, **kwargs) -> None:
        """
        Send an email with the given content to the specified recipient.
        :param to: The email address to send the email to.
        :param subject: The subject of the email.
        :param content_path_html: The relative path to the HTML content of the email where the working directory's default is
                                ``src/casp/services/admin/templates/email``
        :param content_path_plain: The relative path to the plain text content of the email where the working directory's default is
                                ``src/casp/services/admin/templates/email``
        """
        if not self.sg:
            logger.warning("No SendGrid key provided.  Email not sent.")
            return

        with open(f"{self.template_dir}/email/{content_path_html}", "r") as f:
            html_content = f.read().format(**kwargs)
        with open(f"{self.template_dir}/email/{content_path_plain}", "r") as f:
            plain_text_content = f.read().format(**kwargs)

        message = Mail(
            from_email=self.from_address,
            to_emails=to,
            subject=subject,
            html_content=html_content,
            plain_text_content=plain_text_content,
        )
        self.sg.send(message)

    def send_welcome_email(self, email: str, key: str) -> None:
        """
        Send a welcome email to the given email address with the given key.
        :param email: The email address to send the email to.
        :param key: The key to include in the email.
        """
        self._send_email(
            email,
            "Welcome to the Cellarium Cell Annotation Service",
            "welcome.html",
            "welcome.txt",
            key=key,
        )

    def send_new_key_email(self, email: str, key: str) -> None:
        """
        Send an email to the given email address with the given key.
        :param email: The email address to send the email to.
        :param key: The key to include in the email.
        """
        self._send_email(
            email,
            "New Cellarium Cell Annotation Service Key you Requested",
            "new_key.html",
            "new_key.txt",
            key=key,
        )

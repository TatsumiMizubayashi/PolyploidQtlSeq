using System.ComponentModel.DataAnnotations;
using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.QualityControl;
using dq = PolyploidQtlSeqCore.QualityControl;
using op = PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq
{
    /// <summary>
    /// QualityControlコマンド
    /// </summary>
    [Command("qc", Description = "Quality Control by fastp.")]
    internal sealed class QualityControlCommand : CommandBase, IFastpQualityControlCommandOptions
    {
        /// <summary>
        /// QualityControlコマンドを作成する。
        /// </summary>
        public QualityControlCommand()
        {
            InputDir = "";
            OutputDir = "";
            ReadLengthRequired = dq.ReadLengthRequired.DEFAULT;
            NBaseLimit = dq.NBaseLimit.DEFAULT;
            Quality = BaseQuality.DEFAULT;
            CutTailMeanQuality = dq.CutTailMeanQuality.DEFAULT;
            CutTailWindowSize = dq.CutTailWindowSize.DEFAULT;
            ThreadNumber = dq.ThreadNumber.DEFAULT;
            ParameterFile = "";
        }

        [Option(ShortName = InputRawFastqDirectory.SHORT_NAME, LongName = InputRawFastqDirectory.LONG_NAME,
            Description = InputRawFastqDirectory.DESCRIPTION, ValueName = "")]
        public string InputDir { get; set; }

        [Option(ShortName = OutputDirectory.SHORT_NAME, LongName = OutputDirectory.LONG_NAME,
            Description = OutputDirectory.DESCRIPTION, ValueName = "")]
        public string OutputDir { get; set; }

        [Option(ShortName = dq.ReadLengthRequired.SHORT_NAME, LongName = dq.ReadLengthRequired.LONG_NAME,
            Description = dq.ReadLengthRequired.DESCRIPTION, ValueName = "")]
        public int ReadLengthRequired { get; set; }

        [Option(ShortName = dq.NBaseLimit.SHORT_NAME, LongName = dq.NBaseLimit.LONG_NAME,
            Description = dq.NBaseLimit.DESCRIPTION, ValueName = "")]
        public int NBaseLimit { get; set; }

        [Option(ShortName = BaseQuality.SHORT_NAME, LongName = BaseQuality.LONG_NAME,
            Description = BaseQuality.DESCRIPTION, ValueName = "")]
        public int Quality { get; set; }

        [Option(ShortName = dq.CutTailMeanQuality.SHORT_NAME, LongName = dq.CutTailMeanQuality.LONG_NAME,
            Description = dq.CutTailMeanQuality.DESCRIPTION, ValueName = "")]
        public int CutTailMeanQuality { get; set; }

        [Option(ShortName = dq.CutTailWindowSize.SHORT_NAME, LongName = dq.CutTailWindowSize.LONG_NAME, Description = dq.CutTailWindowSize.DESCRIPTION, ValueName = "")]
        public int CutTailWindowSize { get; set; }

        [Option(ShortName = dq.ThreadNumber.SHORT_NAME, LongName = dq.ThreadNumber.LONG_NAME, Description = dq.ThreadNumber.DESCRIPTION, ValueName = "")]
        public int ThreadNumber { get; set; }

        [Option(ShortName = op.ParameterFileParser.SHORT_NAME, LongName = op.ParameterFileParser.LONG_NAME,
            Description = op.ParameterFileParser.DESCRIPTION, ValueName = "")]        
        public string ParameterFile { get; set; }

        public override async Task<int> OnExecuteAsync(CommandLineApplication app)
        {
            int code;
            try
            {
                var fastpQC = new FastpQualityControl(this, app.Options);
                code = await fastpQC.RunAsync();
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine(ex.Message);
                code = 1;
            }

            return code;
        }
    }
}

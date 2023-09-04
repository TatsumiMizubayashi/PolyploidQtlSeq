using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeq.Options.QualityControl;
using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq
{
    /// <summary>
    /// QualityControlコマンド
    /// </summary>
    [Command("qc", Description = "Quality Control by fastp.")]
    internal sealed class QualityControlCommand : CommandBase, IFastpQualityControlSettingValue
    {
        private static readonly string _parameterFileTitle = "Quality Control Parameter file";

        /// <summary>
        /// QualityControlコマンドを作成する。
        /// </summary>
        public QualityControlCommand()
        {
            InputDir = "";
            OutputDir = "";
            ReadLengthRequired = ReadLengthRequiredOption.DEFAULT;
            NBaseLimit = NBaseLimitOption.DEFAULT;
            BaseQuality = BaseQualityOption.DEFAULT;
            CutTailMeanQuality = CutTailMeanQualityOption.DEFAULT;
            CutTailWindowSize = CutTailWindowSizeOption.DEFAULT;
            ThreadNumber = ThreadNumberOption.DEFAULT;
            ParameterFilePath = "";
        }

        [Option(ShortName = InputRawFastqDirectoryOption.SHORT_NAME, LongName = InputRawFastqDirectoryOption.LONG_NAME,
            Description = InputRawFastqDirectoryOption.DESCRIPTION, ValueName = "")]
        public string InputDir { get; set; }

        [Option(ShortName = OutputFastqDirectoryOption.SHORT_NAME, LongName = OutputFastqDirectoryOption.LONG_NAME,
            Description = OutputFastqDirectoryOption.DESCRIPTION, ValueName = "")]
        public string OutputDir { get; set; }

        [Option(ShortName = ReadLengthRequiredOption.SHORT_NAME, LongName = ReadLengthRequiredOption.LONG_NAME,
            Description = ReadLengthRequiredOption.DESCRIPTION, ValueName = "")]
        public int ReadLengthRequired { get; set; }

        [Option(ShortName = NBaseLimitOption.SHORT_NAME, LongName = NBaseLimitOption.LONG_NAME,
            Description = NBaseLimitOption.DESCRIPTION, ValueName = "")]
        public int NBaseLimit { get; set; }

        [Option(ShortName = BaseQualityOption.SHORT_NAME, LongName = BaseQualityOption.LONG_NAME,
            Description = BaseQualityOption.DESCRIPTION, ValueName = "")]
        public int BaseQuality { get; set; }

        [Option(ShortName = CutTailMeanQualityOption.SHORT_NAME, LongName = CutTailMeanQualityOption.LONG_NAME,
            Description = CutTailMeanQualityOption.DESCRIPTION, ValueName = "")]
        public int CutTailMeanQuality { get; set; }

        [Option(ShortName = CutTailWindowSizeOption.SHORT_NAME, LongName = CutTailWindowSizeOption.LONG_NAME, 
            Description = CutTailWindowSizeOption.DESCRIPTION, ValueName = "")]
        public int CutTailWindowSize { get; set; }

        [Option(ShortName = ThreadNumberOption.SHORT_NAME, LongName = ThreadNumberOption.LONG_NAME,
            Description = ThreadNumberOption.DESCRIPTION, ValueName = "")]
        public int ThreadNumber { get; set; }

        [Option(ShortName = ParameterFile.SHORT_NAME, LongName = ParameterFile.LONG_NAME,
            Description = ParameterFile.DESCRIPTION, ValueName = "")]        
        public string ParameterFilePath { get; set; }

        public override async Task<int> OnExecuteAsync(CommandLineApplication app)
        {
            var options = FastpQualityControlOptions.Create(this);
            LoadParameterFile(ParameterFilePath, _parameterFileTitle, options, app);
            if (Validation(options)) Environment.Exit(1);

            int code;
            try
            {
                var fastpQC = new FastpQualityControl(this);
                code = await fastpQC.RunAsync();
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine(ex.Message);
                code = 1;
            }
            finally
            {
                var logParamsFilePath= Path.Combine(OutputDir, "QC.params.txt");
                CreateParameterFile(logParamsFilePath, _parameterFileTitle, options);
            }

            return code;
        }
    }
}

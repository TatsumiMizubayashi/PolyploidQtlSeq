using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.Options;
using PolyploidQtlSeqCore.QualityControl;

namespace PolyploidQtlSeqCore.Application.QualityControl
{
    /// <summary>
    /// Fastp QCコマンドオプション
    /// </summary>
    internal class FastpQualityControlCommandOptions
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [ReadLengthRequired.SHORT_NAME] = ReadLengthRequired.LONG_NAME,
            [ReadLengthRequired.LONG_NAME] = ReadLengthRequired.LONG_NAME,

            [NBaseLimit.SHORT_NAME] = NBaseLimit.LONG_NAME,
            [NBaseLimit.LONG_NAME] = NBaseLimit.LONG_NAME,

            [BaseQuality.SHORT_NAME] = BaseQuality.LONG_NAME,
            [BaseQuality.LONG_NAME] = BaseQuality.LONG_NAME,

            [CutTailMeanQuality.SHORT_NAME] = CutTailMeanQuality.LONG_NAME,
            [CutTailMeanQuality.LONG_NAME] = CutTailMeanQuality.LONG_NAME,

            [CutTailWindowSize.SHORT_NAME] = CutTailWindowSize.LONG_NAME,
            [CutTailWindowSize.LONG_NAME] = CutTailWindowSize.LONG_NAME,

            [ThreadNumber.SHORT_NAME] = ThreadNumber.LONG_NAME,
            [ThreadNumber.LONG_NAME] = ThreadNumber.LONG_NAME,

            [InputRawFastqDirectory.SHORT_NAME] = InputRawFastqDirectory.LONG_NAME,
            [InputRawFastqDirectory.LONG_NAME] = InputRawFastqDirectory.LONG_NAME,

            [OutputDirectory.SHORT_NAME] = OutputDirectory.LONG_NAME,
            [OutputDirectory.LONG_NAME] = OutputDirectory.LONG_NAME
        };

        /// <summary>
        /// Fastp QCコマンドオプションを作成する。
        /// </summary>
        /// <param name="optionValue">オプションの値</param>
        /// <param name="options">オプションリスト</param>
        public FastpQualityControlCommandOptions(IFastpQualityControlCommandOptions optionValue,
            IReadOnlyCollection<CommandOption> options)
        {
            ParameterFile = new ParameterFile(optionValue.ParameterFile);
            var parameterDictionary = ParameterFile.ToParameterDictionary(_toLongNameDictionary);
            var userOptionDictionary = UserSpecifiedLongNameDictionaryCreator.Create(options);

            ReadLengthRequired = new ReadLengthRequired(optionValue.ReadLengthRequired, parameterDictionary, userOptionDictionary);
            NBaseLimit = new NBaseLimit(optionValue.NBaseLimit, parameterDictionary, userOptionDictionary);
            BaseQuality = new BaseQuality(optionValue.Quality, parameterDictionary, userOptionDictionary);
            CutTailMeanQuality = new CutTailMeanQuality(optionValue.CutTailMeanQuality, parameterDictionary, userOptionDictionary);
            CutTailWindowSize = new CutTailWindowSize(optionValue.CutTailWindowSize, parameterDictionary, userOptionDictionary);
            ThreadNumber = new ThreadNumber(optionValue.ThreadNumber, parameterDictionary, userOptionDictionary);
            InputRawFastqDirectory = new InputRawFastqDirectory(optionValue.InputDir, parameterDictionary, userOptionDictionary);
            OutputDirectory = new OutputDirectory(optionValue.OutputDir, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// トリム後リード長の最低値
        /// </summary>
        public ReadLengthRequired ReadLengthRequired { get; }

        /// <summary>
        /// N塩基数の上限値
        /// </summary>
        public NBaseLimit NBaseLimit { get; }

        /// <summary>
        /// 塩基クオリティのしきい値
        /// </summary>
        public BaseQuality BaseQuality { get; }

        /// <summary>
        /// 3'末端トリム時の平均クオリティしきい値
        /// </summary>
        public CutTailMeanQuality CutTailMeanQuality { get; }

        /// <summary>
        /// 3'末端トリム時のウインドウサイズ
        /// </summary>
        public CutTailWindowSize CutTailWindowSize { get; }

        /// <summary>
        /// 使用するスレッド数
        /// </summary>
        public ThreadNumber ThreadNumber { get; }

        /// <summary>
        /// 入力Fastqディレクトリ
        /// </summary>
        public InputRawFastqDirectory InputRawFastqDirectory { get; }

        /// <summary>
        /// 出力ディレクトリ
        /// </summary>
        public OutputDirectory OutputDirectory { get; }

        /// <summary>
        /// パラメータファイル
        /// </summary>
        public ParameterFile ParameterFile { get; }

        /// <summary>
        /// FastpCommonOptionに変換する。
        /// </summary>
        /// <returns>FastpCommonOption</returns>
        public FastpCommonOption ToFastpCommonOption()
        {
            return new FastpCommonOption(
                ReadLengthRequired,
                NBaseLimit,
                BaseQuality,
                CutTailMeanQuality,
                CutTailWindowSize,
                OutputDirectory,
                ThreadNumber);
        }

        /// <summary>
        /// 設定値をパラメータファイル形式で保存する。
        /// </summary>
        /// <param name="filePath">パラメータファイルのPath</param>
        public void SaveParameterFile(string filePath)
        {
            using var writer = new StreamWriter(filePath);

            writer.WriteLine("#Quality Control");
            writer.WriteLine("#LongName\tValue");
            writer.WriteLine(InputRawFastqDirectory.ToParameterFileLine());
            writer.WriteLine(OutputDirectory.ToParameterFileLine());
            writer.WriteLine(ReadLengthRequired.ToParameterFileLine());
            writer.WriteLine(NBaseLimit.ToParameterFileLine());
            writer.WriteLine(BaseQuality.ToParameterFileLine());
            writer.WriteLine(CutTailMeanQuality.ToParameterFileLine());
            writer.WriteLine(CutTailWindowSize.ToParameterFileLine());
            writer.WriteLine(ThreadNumber.ToParameterFileLine());
        }
    }
}

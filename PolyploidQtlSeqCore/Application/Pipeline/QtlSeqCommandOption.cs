using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.Options;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// QtlSeqコマンドオプション
    /// </summary>
    internal class QtlSeqCommandOption
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary;

        /// <summary>
        /// static コンストラクタ
        /// </summary>
        static QtlSeqCommandOption()
        {
            var toLongNameDictionary = new Dictionary<string, string>()
            {
                [ReferenceSequence.SHORT_NAME] = ReferenceSequence.LONG_NAME,
                [ReferenceSequence.LONG_NAME] = ReferenceSequence.LONG_NAME,
            };

            MappingOption.AddLongNameKeyValuePair(toLongNameDictionary);
            AnalysisChrOption.AddLongNameKeyValuePair(toLongNameDictionary);
            BcfToolsVariantCallOption.AddLongNameKeyValuePair(toLongNameDictionary);
            SnpEffOption.AddLongNameKeyValuePair(toLongNameDictionary);
            QtlAnalysisScenarioOptions.AddLongNameKeyValuePair(toLongNameDictionary);

            _toLongNameDictionary = toLongNameDictionary;
        }

        /// <summary>
        /// QtlSeqコマンドオプションを作成する。
        /// </summary>
        /// <param name="optionValues">オプションの値</param>
        /// <param name="options">コマンドオプション</param>
        public QtlSeqCommandOption(IQtlSeqCommandOptions optionValues, IReadOnlyCollection<CommandOption> options)
        {
            ParameterFile = new ParameterFileParser(optionValues.ParameterFile);
            var longNameParameterDictionary = ParameterFile.ToParameterDictionary(_toLongNameDictionary);
            var userOptionDictionary = UserSpecifiedLongNameDictionaryCreator.Create(options);

            ReferenceSequence = new ReferenceSequence(optionValues.ReferenceSequence, longNameParameterDictionary, userOptionDictionary);
            MappingOption = new MappingOption(optionValues, longNameParameterDictionary, userOptionDictionary);
            AnalysisChrOption = new AnalysisChrOption(optionValues, longNameParameterDictionary, userOptionDictionary);
            BcfToolsVariantCallOption = new BcfToolsVariantCallOption(optionValues, longNameParameterDictionary, userOptionDictionary);
            SnpEffOption = new SnpEffOption(optionValues, longNameParameterDictionary, userOptionDictionary);
            QtlAnalysisScenarioOptions = new QtlAnalysisScenarioOptions(optionValues, longNameParameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// リファレンスシークエンスを取得する。
        /// </summary>
        public ReferenceSequence ReferenceSequence { get; }

        /// <summary>
        /// Mappingオプションを取得する。
        /// </summary>
        public MappingOption MappingOption { get; }

        /// <summary>
        /// 解析染色体オプションを取得する。
        /// </summary>
        public AnalysisChrOption AnalysisChrOption { get; }

        /// <summary>
        /// bcftools変異検出オプションを取得する。
        /// </summary>
        public BcfToolsVariantCallOption BcfToolsVariantCallOption { get; }

        /// <summary>
        /// SnpEffオプションを取得する。
        /// </summary>
        public SnpEffOption SnpEffOption { get; }

        /// <summary>
        /// QTL解析シナリオオプションを取得する。
        /// </summary>
        public QtlAnalysisScenarioOptions QtlAnalysisScenarioOptions { get; }

        /// <summary>
        /// パラメータファイルを取得する。
        /// </summary>
        public ParameterFileParser ParameterFile { get; }

        /// <summary>
        /// パラメータファイルを保存する。
        /// </summary>
        /// <param name="filePath">パラメータファイルPath</param>
        public void SaveParameterFile(string filePath)
        {
            using var writer = new StreamWriter(filePath);

            writer.WriteLine("#QTL-Seq Command");
            writer.WriteLine("#LongName\tValue");
            writer.WriteLine(ReferenceSequence.ToParameterFileLine());

            var parameterLineQuery = MappingOption.ToParameterFileLines()
                .Concat(AnalysisChrOption.ToParameterFileLines())
                .Concat(BcfToolsVariantCallOption.ToParameterFileLines())
                .Concat(SnpEffOption.ToParameterFileLines())
                .Concat(QtlAnalysisScenarioOptions.ToParameterFileLines());
            foreach (var line in parameterLineQuery)
            {
                writer.WriteLine(line);
            }
        }
    }
}

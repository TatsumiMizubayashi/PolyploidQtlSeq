using McMaster.Extensions.CommandLineUtils;

namespace PolyploidQtlSeqCore.Options
{
    /// <summary>
    /// オプションコレクションの抽象クラス
    /// </summary>
    public abstract class OptionCollection
    {
        private readonly Option[] _options;

        /// <summary>
        /// オプションコレクション インスタンスを作成する。
        /// </summary>
        /// <param name="options">オプション配列</param>
        public OptionCollection(Option[] options)
        {
            _options = options;
        }

        /// <summary>
        /// オプション名に該当するLongNameを取得する。
        /// 該当するLongNameがない場合は空文字を返す。
        /// </summary>
        /// <param name="name">オプション名</param>
        /// <returns>LongName</returns>
        public string GetLongName(string name)
        {
            var longNames = _options.Select(x => x.GetLongName(name))
                .Where(x => !string.IsNullOrEmpty(x))
                .ToArray();

            return longNames.Length == 1 ? longNames[0] : "";
        }

        /// <summary>
        /// オプション値を設定する。
        /// ユーザー指定オプションが指定されている場合はparamTableの値は設定しない。
        /// </summary>
        /// <param name="longNameParamsDictionary">LongNameパラメーター辞書</param>
        /// <param name="commandOptions">ユーザー指定オプション</param>
        public void SetValues(IReadOnlyDictionary<string, string> longNameParamsDictionary, IReadOnlyCollection<CommandOption> commandOptions)
        {
            var userSetLongNameTable = UserSpecifiedLongNameDictionaryCreator.Create(commandOptions);

            foreach (var option in _options)
            {
                option.SetValue(longNameParamsDictionary, userSetLongNameTable);
            }
        }

        /// <summary>
        /// データ検証結果を取得する。
        /// </summary>
        /// <returns>エラーがあるデータ検証結果</returns>
        public DataValidationResult[] Validation()
        {
            var varidationResults = _options.Select(x => x.Validation())
                .Where(x => x.HasError)
                .ToArray();

            return varidationResults;
        }

        /// <summary>
        /// Key-Value値を取得する。
        /// </summary>
        /// <returns>Key-Value値</returns>
        public string[] GetKeyValues()
        {
            return [.. _options.Select(x => x.GetKeyValue())];
        }
    }
}

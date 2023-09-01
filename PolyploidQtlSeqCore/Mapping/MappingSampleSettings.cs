namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Mappingサンプル設定
    /// </summary>
    [Obsolete("オプションスイッチ機能を削除予定")]
    internal class MappingSampleSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [Parent1Directory.SHORT_NAME] = Parent1Directory.LONG_NAME,
            [Parent1Directory.LONG_NAME] = Parent1Directory.LONG_NAME,

            [Parent2Directory.SHORT_NAME] = Parent2Directory.LONG_NAME,
            [Parent2Directory.LONG_NAME] = Parent2Directory.LONG_NAME,

            [Bulk1Directory.SHORT_NAME] = Bulk1Directory.LONG_NAME,
            [Bulk1Directory.LONG_NAME] = Bulk1Directory.LONG_NAME,

            [Bulk2Directory.SHORT_NAME] = Bulk2Directory.LONG_NAME,
            [Bulk2Directory.LONG_NAME] = Bulk2Directory.LONG_NAME
        };

        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }

        /// <summary>
        /// Mappngオプションを作成する。
        /// </summary>
        /// <param name="settingValue">Mappingサンプル設定値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定オプション辞書</param>
        public MappingSampleSettings(IMappingSampleSettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Parent1Directory = new Parent1Directory(settingValue.Parent1Dir, parameterDictionary, userOptionDictionary);
            Parent2Directory = new Parent2Directory(settingValue.Parent2Dir, parameterDictionary, userOptionDictionary);
            Bulk1Directory = new Bulk1Directory(settingValue.Bulk1Dir, parameterDictionary, userOptionDictionary);
            Bulk2Directory = new Bulk2Directory(settingValue.Bulk2Dir, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// Parent1ディレクトリ
        /// </summary>
        public Parent1Directory Parent1Directory { get; }

        /// <summary>
        /// Parent2ディレクトリ
        /// </summary>
        public Parent2Directory Parent2Directory { get; }

        /// <summary>
        /// Bulk1ディレクトリ
        /// </summary>
        public Bulk1Directory Bulk1Directory { get; }

        /// <summary>
        /// Bulk2ディレクトリ
        /// </summary>
        public Bulk2Directory Bulk2Directory { get; }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                Parent1Directory.ToParameterFileLine(),
                Parent2Directory.ToParameterFileLine(),
                Bulk1Directory.ToParameterFileLine(),
                Bulk2Directory.ToParameterFileLine()
            };
        }
    }
}
